# inputs: development polygon, population raster
#

import os
import arcpy

# Enable GDAL memory datasets (required for rasterio.mask in newer GDAL versions)
os.environ["GDAL_MEM_ENABLE_OPEN"] = "YES"

import pyproj
from shapely import Polygon
from shapely.geometry import box, shape
from utils.api import fetch_isochrone, get_parks_gdf_via_arc
from constants import DEVELOPMENT_POLYGONS_GEOJSON
import geopandas as gpd
from tqdm import tqdm
import rasterio
import numpy as np
from contextlib import nullcontext
from utils.arc_utils import gdf_to_featureclass
from utils.utils import (
    calculate_demand_metrics_2sfca,
    distribute_population_stats,
    get_raster_clip_under_polygon,
    overlay_rasters,
    create_memory_raster,
)

DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"

development_polygon_one = Polygon(
    DEVELOPMENT_POLYGONS_GEOJSON["features"][0]["geometry"]["coordinates"][0]
)
development_polygon_two = Polygon(
    DEVELOPMENT_POLYGONS_GEOJSON["features"][1]["geometry"]["coordinates"][0]
)


def append_isochrones_to_parks_gdf(parks_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    isochrones = []
    # Convert parks to WGS84 for API calls
    parks_wgs84 = parks_gdf.to_crs("EPSG:4326")
    for _, row in tqdm(
        parks_wgs84.iterrows(),
        total=len(parks_wgs84),
        desc="Fetching isochrones for parks",
    ):
        centroid = row["geometry"].centroid
        isochrone_geojson = fetch_isochrone(centroid.x, centroid.y)

        if isochrone_geojson:
            isochrone_polygon = shape(isochrone_geojson["geometry"])
            isochrones.append(isochrone_polygon)
        else:
            isochrones.append(None)

    parks_gdf["isochrone_polygon"] = gpd.GeoSeries(isochrones, crs=parks_wgs84.crs)
    return parks_gdf


def append_demand_metrics(
    parks_gdf: gpd.GeoDataFrame,
    population_raster: str | rasterio.DatasetReader,
    demand_field_prefix: str,
    decay_type: str = "gaussian",
) -> gpd.GeoDataFrame:
    """Calculate 2SFCA demand metrics for each park's isochrone.

    Args:
        parks_gdf: GeoDataFrame with park geometries and isochrone_polygon column
        population_raster: Either a file path string or an open rasterio DatasetReader
        demand_field_prefix: Prefix for output columns (e.g., "current" -> "current_pop", "current_m2_per_person", etc.)
        decay_type: Distance decay function ("none", "inverse", "inverse_square", "gaussian")

    Output columns added:
        {prefix}_pop: Raw population sum in isochrone
        {prefix}_pop_weighted: Distance-weighted population sum
        {prefix}_park_area_m2: Park area in square meters
        {prefix}_m2_per_person: Square meters of park per person (2SFCA ratio)
        {prefix}_m2_per_person_weighted: Same but with distance decay
        {prefix}_acres_per_1000: Acres per 1000 people
        {prefix}_acres_per_1000_weighted: Same but with distance decay
    """
    # Handle both file path and already-open raster
    ctx = (
        rasterio.open(population_raster)
        if isinstance(population_raster, str)
        else nullcontext(population_raster)
    )

    with ctx as src:
        # Transformer for isochrone polygons (WGS84 -> raster CRS)
        isochrone_transformer = pyproj.Transformer.from_crs(
            parks_gdf[
                "isochrone_polygon"
            ].crs,  # Source CRS from isochrone polygons (WGS84)
            src.crs,  # Target CRS from raster
            always_xy=True,
        )

        # Transformer for park geometries (Web Mercator -> raster CRS)
        park_transformer = pyproj.Transformer.from_crs(
            parks_gdf.crs,  # Source CRS from park geometries (typically Web Mercator)
            src.crs,  # Target CRS from raster
            always_xy=True,
        )

        for idx, row in tqdm(
            parks_gdf.iterrows(),
            total=len(parks_gdf),
            desc=f"Calculating {demand_field_prefix} demand metrics",
        ):
            metrics = calculate_demand_metrics_2sfca(
                isochrone_polygon=row["isochrone_polygon"],
                park_geometry=row["geometry"],
                raster=src,
                isochrone_transformer=isochrone_transformer,
                park_transformer=park_transformer,
                decay_type=decay_type,
            )

            # Add each metric as a column with the prefix
            for metric_name, value in metrics.items():
                col_name = f"{demand_field_prefix}_{metric_name}"
                parks_gdf.at[idx, col_name] = value

    return parks_gdf


def demand_analysis(developments_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Analyze park demand impact from multiple development polygons.

    Args:
        developments_gdf: GeoDataFrame with columns:
            - geometry: Development polygon geometries (should be in WGS84)
            - population_change: Float column with population change for each polygon

    Returns:
        GeoDataFrame of parks with current and projected demand metrics
    """
    # Validate input
    if "population_change" not in developments_gdf.columns:
        raise ValueError("developments_gdf must have a 'population_change' column")

    # Ensure developments are in WGS84
    if developments_gdf.crs is None:
        print("Warning: developments_gdf has no CRS, assuming WGS84")
        developments_gdf = developments_gdf.set_crs("EPSG:4326")
    elif developments_gdf.crs != "EPSG:4326":
        developments_gdf = developments_gdf.to_crs("EPSG:4326")

    # Get isochrones for all development polygons to determine search area
    print(f"Processing {len(developments_gdf)} development polygons...")

    development_isochrones = []
    for idx, row in tqdm(
        developments_gdf.iterrows(),
        total=len(developments_gdf),
        desc="Fetching development isochrones",
    ):
        centroid = row.geometry.centroid
        isochrone_geojson = fetch_isochrone(centroid.x, centroid.y)
        if isochrone_geojson:
            development_isochrones.append(shape(isochrone_geojson["geometry"]))
        else:
            # Fallback: use a buffer around the polygon
            development_isochrones.append(
                row.geometry.buffer(0.01)
            )  # ~1km buffer in degrees

    developments_gdf["isochrone"] = gpd.GeoSeries(
        development_isochrones, crs="EPSG:4326"
    )

    # Get combined bounds of all development isochrones for park search
    from shapely.ops import unary_union

    combined_isochrone = unary_union(development_isochrones)
    combined_centroid = combined_isochrone.centroid

    # Calculate search radius from combined isochrone bounds
    radius_degrees = combined_centroid.hausdorff_distance(combined_isochrone.boundary)
    radius_meters = radius_degrees * 111000
    print(
        f"Searching for parks within {radius_meters:.0f} meters of development area center"
    )

    # Find parks in the combined area
    parks_gdf = get_parks_gdf_via_arc(
        combined_centroid, radius_meters, point_crs="EPSG:4326"
    )

    if len(parks_gdf) == 0:
        print("No parks found in the search area")
        return parks_gdf

    parks_gdf = append_isochrones_to_parks_gdf(parks_gdf)
    parks_gdf = append_demand_metrics(parks_gdf, DEFAULT_GHS_RASTER, "current")

    # Set up the projected demand raster
    all_isochrone_bounds = parks_gdf["isochrone_polygon"].total_bounds
    all_isochrone_bounds_polygon = box(*all_isochrone_bounds)

    isochrone_crs = parks_gdf["isochrone_polygon"].crs

    with rasterio.open(DEFAULT_GHS_RASTER) as population_raster:
        transformer = pyproj.Transformer.from_crs(
            isochrone_crs, population_raster.crs, always_xy=True
        )
        raster_crs = population_raster.crs
        raster_nodata = population_raster.nodata

        # Clip population raster to study area (base layer)
        population_array, _, population_transform = get_raster_clip_under_polygon(
            all_isochrone_bounds_polygon, population_raster, transformer
        )

        # Start with the base population
        projected_population_array = population_array.copy().astype(np.float32)

        # Overlay each development polygon's population onto the base
        for idx, row in tqdm(
            developments_gdf.iterrows(),
            total=len(developments_gdf),
            desc="Distributing development population",
        ):
            development_polygon = row.geometry
            population_change = row["population_change"]

            if population_change == 0:
                continue

            # Distribute population for this development
            development_array, development_transform = distribute_population_stats(
                development_polygon, population_change, population_raster, transformer
            )

            # Overlay onto the projected population
            projected_population_array = overlay_rasters(
                projected_population_array,
                population_transform,
                development_array,
                development_transform,
            )

    # Calculate projected demand using the combined population raster
    with create_memory_raster(
        projected_population_array, population_transform, raster_crs, raster_nodata
    ) as projected_raster:
        parks_gdf = append_demand_metrics(parks_gdf, projected_raster, "projected")

    # Calculate change from current to projected (positive = increase)
    parks_gdf["pop_change"] = parks_gdf["projected_pop"] - parks_gdf["current_pop"]
    parks_gdf["m2_per_person_change"] = (
        parks_gdf["projected_m2_per_person"] - parks_gdf["current_m2_per_person"]
    )
    parks_gdf["acres_per_1000_change"] = (
        parks_gdf["projected_acres_per_1000"] - parks_gdf["current_acres_per_1000"]
    )

    # Percent change (useful for comparing across different sized parks)
    parks_gdf["pop_pct_change"] = (
        parks_gdf["pop_change"] / parks_gdf["current_pop"]
    ) * 100
    parks_gdf["acres_per_1000_pct_change"] = (
        parks_gdf["acres_per_1000_change"] / parks_gdf["current_acres_per_1000"]
    ) * 100

    return parks_gdf


# Create a GeoDataFrame with development polygons and their population changes
developments_gdf = gpd.GeoDataFrame(
    {
        "geometry": [development_polygon_one, development_polygon_two],
        "population_change": [5713.0, 4205.0],
    },
    crs="EPSG:4326",
)

# Run the analysis
arcpy.AddMessage("Starting demand analysis...")
parks_gdf = demand_analysis(developments_gdf)

# Add results to map
arcpy.AddMessage(f"Analysis complete. Adding {len(parks_gdf)} parks to map...")

# Save to a shapefile instead of in-memory (more reliable)
output_shp = r"C:\gispy\scratch\parks_demand_analysis.shp"

fc_path = gdf_to_featureclass(parks_gdf, output_path=output_shp)
arcpy.AddMessage(f"Saved to: {fc_path}")

# Add to map
try:
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    m = aprx.activeMap
    if m:
        arcpy.AddMessage(f"Active map: {m.name}")
        m.addDataFromPath(fc_path)
        arcpy.AddMessage("Layer added!")
    else:
        arcpy.AddWarning("No active map found")
except Exception as e:
    arcpy.AddError(f"Error adding to map: {e}")
