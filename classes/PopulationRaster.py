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
        self,
        polygon: Polygon,
        polygon_crs=None,
        all_touched: bool = True
    ) -> tuple[np.ndarray, np.ndarray, Affine]:
        """
        Clip the raster to a polygon.
        
        Args:
            polygon: Shapely polygon geometry
            polygon_crs: CRS of the polygon (e.g., "EPSG:4326"). If None, assumes raster CRS.
            all_touched: If True, include all pixels touched by polygon
            
        Returns:
            Tuple of (clipped_array, valid_mask, transform) where:
            - clipped_array: The clipped raster data
            - valid_mask: Boolean mask indicating non-nodata pixels
            - transform: Affine transform for the clipped raster
        """
        # Transform polygon to raster CRS if needed
        if polygon_crs is not None:
            transformer = self.get_transformer(polygon_crs)
            polygon = shapely_transform(transformer.transform, polygon)
        
        out_image, out_transform = mask(
            self._src, [polygon], crop=True, all_touched=all_touched
        )
        
        raster_array = out_image[0]
        valid_mask = raster_array != self.nodata
        
        return raster_array, valid_mask, out_transform
    
    def distribute_population(
        self,
        polygon: Polygon,
        population_change: float,
        polygon_crs=None
    ) -> tuple[np.ndarray, Affine]:
        """
        Distribute population change randomly across valid cells in a polygon.
        
        Args:
            polygon: Shapely polygon (e.g., a development area)
            population_change: Total population to distribute
            polygon_crs: CRS of the polygon. If None, assumes raster CRS.
            
        Returns:
            Tuple of (output_array, transform) containing the distributed population
        """
        raster_array, valid_mask, out_transform = self.clip_to_polygon(
            polygon, polygon_crs
        )
        
        valid_cells = np.where(valid_mask)
        num_valid_cells = len(valid_cells[0])
        
        if num_valid_cells == 0:
            raise ValueError("No valid cells found in the polygon")
        
        # Distribute population randomly across valid cells
        random_weights = np.random.rand(num_valid_cells)
        random_weights /= np.sum(random_weights)
        
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
        sigma_km: float = 0.5
    ) -> dict:
        """
        Calculate 2SFCA demand metrics for a park.
        
        Args:
            isochrone_polygon: The catchment area (walking zone) for the park
            park_geometry: The park polygon (for area and centroid)
            isochrone_crs: CRS of the isochrone. If None, assumes raster CRS.
            park_crs: CRS of the park geometry. If None, assumes raster CRS.
            decay_type: "none", "inverse", "inverse_square", or "gaussian"
            sigma_km: For gaussian decay, controls spread (higher = slower decay)
            
        Returns:
            Dictionary with metrics:
            - pop: Raw population sum in isochrone
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
                park_in_raster_crs = shapely_transform(transformer.transform, park_geometry)
            else:
                park_in_raster_crs = park_geometry
            
            # Apply distance decay
            weighted_pop = self._apply_distance_decay(
                raster_array, out_transform,
                park_centroid.x, park_centroid.y,
                decay_type, sigma_km
            )
            
            # Calculate metrics
            pop_raw = float(np.sum(raster_array[valid_mask]))
            pop_weighted = float(np.sum(weighted_pop[valid_mask]))
            park_area_m2 = park_in_raster_crs.area
            
            # 2SFCA ratios
            m2_per_person = park_area_m2 / pop_raw if pop_raw > 0 else None
            m2_per_person_weighted = park_area_m2 / pop_weighted if pop_weighted > 0 else None
            
            # Acres per 1000 people (more readable)
            acres_per_1000 = (park_area_m2 / 4047) / (pop_raw / 1000) if pop_raw > 0 else None
            acres_per_1000_weighted = (park_area_m2 / 4047) / (pop_weighted / 1000) if pop_weighted > 0 else None
            
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
        sigma_km: float = 0.5
    ) -> np.ndarray:
        """Weight population cells by distance to a target point."""
        if decay_type == "none":
            return pop_array
        
        rows, cols = pop_array.shape
        col_indices, row_indices = np.meshgrid(np.arange(cols), np.arange(rows))
        
        # Convert pixel indices to world coordinates (cell centers)
        cell_x = transform.c + (col_indices + 0.5) * transform.a
        cell_y = transform.f + (row_indices + 0.5) * transform.e
        
        # Calculate distances in km
        distances_km = np.sqrt((cell_x - target_x)**2 + (cell_y - target_y)**2) / 1000
        
        # Apply decay function
        if decay_type == "inverse":
            weights = 1 / (1 + distances_km)
        elif decay_type == "inverse_square":
            weights = 1 / (1 + distances_km**2)
        elif decay_type == "gaussian":
            weights = np.exp(-(distances_km**2) / (2 * sigma_km**2))
        else:
            weights = np.ones_like(distances_km)
        
        return pop_array * weights
    
    def __enter__(self):
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._owns_file:
            self._src.close()
        return False