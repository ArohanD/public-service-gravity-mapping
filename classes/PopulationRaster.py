import rasterio
from rasterio import Affine
from rasterio.mask import mask
import numpy as np
from shapely import Polygon
from shapely.ops import transform as shapely_transform
import pyproj


class PopulationRaster:
    """
    Wraps a population raster with coordinate transformation and analysis utilities.

    Handles CRS transformations internally, so you just pass geometries in their
    native CRS and the class handles reprojection to the raster's CRS.

    ## Two-Step Floating Catchment Area (2SFCA) Method

    This class implements 2SFCA for measuring park accessibility. The key idea:

    **Step 1 (Supply Ratio):** For each park, calculate how much space is available
    per person who can reach it:

        R = Park Area / Population in Catchment (isochrone for park)

    **Step 2 (Accessibility):** For each location (pixel), sum up the ratios of all parks
    they can reach to get an accessibility score.

    ## How Population is Counted

    We use a population raster where each pixel (e.g., 100m x 100m) contains a population
    count. The catchment is defined by an isochrone polygon (everywhere reachable within
    X minutes).

        Population Raster           Clipped to Isochrone
        ┌───┬───┬───┬───┐          ┌───┬───┬───┬───┐
        │ 5 │12 │ 8 │ 2 │          │   │12 │ 8 │   │
        ├───┼───┼───┼───┤   -->    ├───┼───┼───┼───┤
        │20 │35 │15 │ 4 │          │20 │35 │15 │   │
        └───┴───┴───┴───┘          └───┴───┴───┴───┘

        Total population in catchment = 12 + 8 + 20 + 35 + 15 = 90

    ## Distance Decay Weighting

    Basic 2SFCA treats everyone equally. But someone 1 minute away uses the park
    more than someone 14 minutes away. Distance decay fixes this by weighting
    each pixel's population by how far it is from the park:

        weighted_pop = population × weight(distance)

    Available decay functions are listed in _apply_distance_decay.

    Usage:
        with PopulationRaster("path/to/raster.tif") as pop:
            # Clip to a polygon (in WGS84)
            array, valid_mask, transform = pop.clip_to_polygon(polygon, "EPSG:4326")

            # Calculate demand metrics for a park
            metrics = pop.calculate_demand_metrics(
                isochrone_polygon=isochrone,
                park_geometry=park,
                isochrone_crs="EPSG:4326",
                park_crs="EPSG:3857"
            )

            # Distribute population change to a development polygon
            array, transform = pop.distribute_population(polygon, 500, "EPSG:4326")
    """

    def __init__(self, source):
        """
        Initialize from a file path or existing rasterio DatasetReader.

        Args:
            source: Either a file path string or an open rasterio.DatasetReader
        """
        if isinstance(source, str):
            self._src = rasterio.open(source)
            self._owns_file = True
        else:
            self._src = source
            self._owns_file = False

        # Cache for transformers (source_crs -> transformer)
        self._transformers = {}

    @property
    def crs(self):
        """The raster's coordinate reference system."""
        return self._src.crs

    @property
    def nodata(self):
        """The raster's nodata value."""
        return self._src.nodata

    @property
    def transform(self):
        """The raster's affine transform."""
        return self._src.transform

    def get_transformer(self, source_crs) -> pyproj.Transformer:
        """
        Get a transformer FROM source_crs TO this raster's CRS.

        Transformers are cached for reuse.
        """
        crs_key = str(source_crs)
        if crs_key not in self._transformers:
            self._transformers[crs_key] = pyproj.Transformer.from_crs(
                source_crs, self.crs, always_xy=True
            )
        return self._transformers[crs_key]

    def clip_to_polygon(
        self, polygon: Polygon, polygon_crs=None, all_touched: bool = True
    ) -> tuple[np.ndarray, np.ndarray, Affine]:
        """
        Clip the raster to a polygon. This is used for both the isochrone and
        the development polygons.

        In 2SFCA, this is how we define "who can reach the park." The polygon
        is typically an isochrone (everywhere reachable within X minutes), and
        clipping the population raster to it gives us only the relevant pixels.

        Only pixels inside the isochrone are kept and the raster is resized accordingly.

        Args:
            polygon: Shapely polygon geometry
            polygon_crs: CRS of the polygon (e.g., "EPSG:4326"). If None, assumes raster CRS.
            all_touched: If True, include all pixels touched by polygon edge
                         (recommended for catchment analysis to avoid gaps)

        Returns:
            Tuple of (clipped_array, valid_mask, transform) where:
            - clipped_array: The clipped raster data (population values)
            - valid_mask: Boolean mask indicating non-nodata pixels (where people live)
            - transform: Affine transform for the clipped raster (for coordinate math)
        """
        # Transform polygon to raster CRS if needed
        # (polygon might be in lat/lon, raster might be in meters)
        if polygon_crs is not None:
            transformer = self.get_transformer(polygon_crs)
            polygon = shapely_transform(transformer.transform, polygon)

        # Use rasterio's mask function to clip and crop in one step
        out_image, out_transform = mask(
            self._src, [polygon], crop=True, all_touched=all_touched
        )

        raster_array = out_image[0]
        # valid_mask tells us which pixels have actual population data
        # (vs nodata pixels like water, outside coverage, etc.)
        valid_mask = raster_array != self.nodata

        return raster_array, valid_mask, out_transform

    def distribute_population(
        self, polygon: Polygon, population_change: float, polygon_crs=None
    ) -> tuple[np.ndarray, Affine]:
        """
        Distribute population change across valid cells in a polygon.

        This is used to model future development: given a development polygon
        and the expected number of new residents, distribute that population
        across the raster cells within the polygon.

        Example: A new apartment complex (polygon) will add 500 people.

            Development Polygon          Distributed Population
            ┌───────────────┐            ┌───┬───┬───┐
            │               │            │ 85│120│ 95│  ← Random distribution
            │   500 people  │    -->     ├───┼───┼───┤    totaling 500
            │               │            │110│ 90│   │
            └───────────────┘            └───┴───┴───┘

        The output array can be overlaid onto the base population raster to
        create a "projected" population raster for future demand analysis.

        Args:
            polygon: Shapely polygon (e.g., a development area)
            population_change: Total population to distribute
            polygon_crs: CRS of the polygon. If None, assumes raster CRS.

        Returns:
            Tuple of (output_array, transform) where:
            - output_array: Raster with distributed population (same grid as input)
            - transform: Affine transform for positioning this array
        """
        # Clip to get the valid cells within the development polygon
        raster_array, valid_mask, out_transform = self.clip_to_polygon(
            polygon, polygon_crs
        )

        valid_cells = np.where(valid_mask)
        num_valid_cells = len(valid_cells[0])

        if num_valid_cells == 0:
            raise ValueError("No valid cells found in the polygon")

        # Distribute population randomly across valid cells
        # Each cell gets a random weight, then we normalize so weights sum to 1
        random_weights = np.random.rand(num_valid_cells)
        random_weights /= np.sum(random_weights)

        # Create output array and assign population proportionally
        output = np.zeros(raster_array.shape, dtype=np.float32)
        output[valid_cells] = random_weights * population_change

        return output, out_transform

    def calculate_demand_metrics(
        self,
        isochrone_polygon: Polygon,
        park_geometry: Polygon,
        isochrone_crs=None,
        park_crs=None,
        decay_type: str = "gaussian",
        sigma_km: float = 0.5,
    ) -> dict:
        """
        Calculate 2SFCA demand metrics for a park.

        This implements Step 1 of the Two-Step Floating Catchment Area method:
        calculate the supply-to-demand ratio for a single park.

            R = Park Area / Population in Catchment

        With distance decay, population is weighted so closer people count more:

            R_weighted = Park Area / Σ(population × distance_weight)

        Example: A 5-acre park with 1000 people in its catchment:
            - acres_per_1000 = 5 / (1000/1000) = 5.0 acres per 1000 people
            - If weighted pop = 600, acres_per_1000_weighted = 5 / 0.6 = 8.3
              (higher because distant people count less)

        Args:
            isochrone_polygon: The catchment area (travel time zone) for the park.
                               This defines who can "reach" the park.
            park_geometry: The park polygon (for area and centroid)
            isochrone_crs: CRS of the isochrone. If None, assumes raster CRS.
            park_crs: CRS of the park geometry. If None, assumes raster CRS.
            decay_type: "none", "inverse", "inverse_square", or "gaussian"
            sigma_km: For gaussian decay, controls spread (higher = slower decay)

        Returns:
            Dictionary with metrics:
            - pop: Raw population sum in isochrone (unweighted)
            - pop_weighted: Distance-weighted population sum
            - park_area_m2: Park area in square meters
            - m2_per_person: Square meters per person
            - m2_per_person_weighted: Same with distance decay
            - acres_per_1000: Acres per 1000 people
            - acres_per_1000_weighted: Same with distance decay
        """
        try:
            # Clip population to isochrone
            raster_array, valid_mask, out_transform = self.clip_to_polygon(
                isochrone_polygon, isochrone_crs, all_touched=True
            )

            # Transform park centroid to raster CRS for distance calculations
            park_centroid = park_geometry.centroid
            if park_crs is not None:
                transformer = self.get_transformer(park_crs)
                park_centroid = shapely_transform(transformer.transform, park_centroid)
                park_in_raster_crs = shapely_transform(
                    transformer.transform, park_geometry
                )
            else:
                park_in_raster_crs = park_geometry

            # Apply distance decay
            weighted_pop = self._apply_distance_decay(
                raster_array,
                out_transform,
                park_centroid.x,
                park_centroid.y,
                decay_type,
                sigma_km,
            )

            # Calculate metrics
            pop_raw = float(np.sum(raster_array[valid_mask]))
            pop_weighted = float(np.sum(weighted_pop[valid_mask]))
            park_area_m2 = park_in_raster_crs.area

            # 2SFCA supply ratio: R = Supply / Demand
            # Higher values = more space per person = less crowded
            m2_per_person = park_area_m2 / pop_raw if pop_raw > 0 else None
            m2_per_person_weighted = (
                park_area_m2 / pop_weighted if pop_weighted > 0 else None
            )

            # Conversion: 1 acre = 4047 m²
            acres_per_1000 = (
                (park_area_m2 / 4047) / (pop_raw / 1000) if pop_raw > 0 else None
            )
            acres_per_1000_weighted = (
                (park_area_m2 / 4047) / (pop_weighted / 1000)
                if pop_weighted > 0
                else None
            )

            return {
                "pop": pop_raw,
                "pop_weighted": pop_weighted,
                "park_area_m2": park_area_m2,
                "m2_per_person": m2_per_person,
                "m2_per_person_weighted": m2_per_person_weighted,
                "acres_per_1000": acres_per_1000,
                "acres_per_1000_weighted": acres_per_1000_weighted,
            }

        except Exception as e:
            print(f"Error calculating demand metrics: {e}")
            return {
                "pop": None,
                "pop_weighted": None,
                "park_area_m2": None,
                "m2_per_person": None,
                "m2_per_person_weighted": None,
                "acres_per_1000": None,
                "acres_per_1000_weighted": None,
            }

    def _apply_distance_decay(
        self,
        pop_array: np.ndarray,
        transform: Affine,
        target_x: float,
        target_y: float,
        decay_type: str = "gaussian",
        sigma_km: float = 0.5,
    ) -> np.ndarray:
        """
        Weight population cells by distance to a target point (e.g., park center).

        Distance decay addresses a limitation of basic 2SFCA: it treats everyone
        in the catchment equally, but in reality, closer people use a park more.

            Basic 2SFCA:    Person 1 min away  = Person 14 min away  (both count as 1)
            With decay:     Person 1 min away >> Person 14 min away  (weighted)

        Each pixel's population is multiplied by a weight between 0 and 1:
            - weight ≈ 1.0 → very close, full count
            - weight ≈ 0.5 → medium distance, half count
            - weight ≈ 0.0 → far away, barely counts

        Available decay functions and their characteristics:

            Distance (km)  │ Inverse │ Inverse² │ Gaussian (σ=0.5)
            ───────────────┼─────────┼──────────┼──────────────────
                 0.0       │  1.00   │   1.00   │   1.00
                 0.25      │  0.80   │   0.94   │   0.88
                 0.5       │  0.67   │   0.80   │   0.61
                 1.0       │  0.50   │   0.50   │   0.14
                 2.0       │  0.33   │   0.20   │   0.00

        - Gaussian: Steep drop-off, emphasizes nearby population (recommended)
        - Inverse: Gentle decay, distant people still count somewhat
        - Inverse Square: Middle ground between the two

        Args:
            pop_array: 2D array of population values (one per pixel)
            transform: Affine transform to convert pixel indices to coordinates
            target_x, target_y: Coordinates of the target point (park center)
            decay_type: "none", "inverse", "inverse_square", or "gaussian"
            sigma_km: For gaussian, controls the spread (higher = slower decay)
                      σ = 0.5 means weight drops to ~61% at 0.5 km

        Returns:
            Array of same shape with population × weight for each pixel
        """
        if decay_type == "none":
            return pop_array

        rows, cols = pop_array.shape
        col_indices, row_indices = np.meshgrid(np.arange(cols), np.arange(rows))

        # Convert pixel indices to world coordinates (cell centers)
        # The +0.5 shifts from corner to center of each pixel
        cell_x = transform.c + (col_indices + 0.5) * transform.a
        cell_y = transform.f + (row_indices + 0.5) * transform.e

        # Calculate Euclidean distance from each pixel center to park center
        # Euclidian distance is not the most realistic distance metric, but the
        # density of pixels and the emphasis on closer pixels being weighted more
        # makes it a reasonable approximation.
        # Divide by 1000 to convert meters to kilometers
        distances_km = (
            np.sqrt((cell_x - target_x) ** 2 + (cell_y - target_y) ** 2) / 1000
        )

        # Apply decay function
        if decay_type == "inverse":
            # weight = 1 / (1 + d)
            # Gentle decay: at 1 km, weight = 0.5
            weights = 1 / (1 + distances_km)
        elif decay_type == "inverse_square":
            # weight = 1 / (1 + d²)
            # Steeper than inverse: at 1 km, weight = 0.5
            weights = 1 / (1 + distances_km**2)
        elif decay_type == "gaussian":
            # weight = e^(-d² / 2σ²)
            # Bell curve: at d = σ, weight ≈ 0.61
            # This is the most commonly used in 2SFCA literature
            weights = np.exp(-(distances_km**2) / (2 * sigma_km**2))
        else:
            weights = np.ones_like(distances_km)

        return pop_array * weights

    # Following methods allow us to use "with" syntax to open and close the raster.
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._owns_file:
            self._src.close()
        return False
