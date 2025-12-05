# inputs: development polygon, population raster
#

import os
# Enable GDAL memory datasets (required for rasterio.mask in newer GDAL versions)
os.environ['GDAL_MEM_ENABLE_OPEN'] = 'YES'

import pyproj
from shapely import Polygon
from shapely.geometry import box, shape
from get_isochrones import fetch_isochrone
from get_parks import get_parks_gdf_via_arc
from constants import DEVELOPMENT_POLYGON_GEOJSON
import geopandas as gpd
from tqdm import tqdm
import rasterio
import numpy as np
from contextlib import nullcontext

from utils.utils import calculate_demand_metrics_2sfca, distribute_population_stats, get_raster_clip_under_polygon, overlay_rasters, create_memory_raster

DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"

development_polygon = Polygon(DEVELOPMENT_POLYGON_GEOJSON["features"][0]["geometry"]["coordinates"][0])

def append_isochrones_to_parks_gdf(parks_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    isochrones = []
    # Convert parks to WGS84 for API calls
    parks_wgs84 = parks_gdf.to_crs("EPSG:4326")
    for _, row in tqdm(parks_wgs84.iterrows(), total=len(parks_wgs84), desc="Fetching isochrones for parks"):
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
    decay_type: str = "gaussian"
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
    ctx = rasterio.open(population_raster) if isinstance(population_raster, str) else nullcontext(population_raster)
    
    with ctx as src:
        # Transformer for isochrone polygons (WGS84 -> raster CRS)
        isochrone_transformer = pyproj.Transformer.from_crs(
            parks_gdf["isochrone_polygon"].crs,  # Source CRS from isochrone polygons (WGS84)
            src.crs,  # Target CRS from raster
            always_xy=True
        )
        
        # Transformer for park geometries (Web Mercator -> raster CRS)
        park_transformer = pyproj.Transformer.from_crs(
            parks_gdf.crs,  # Source CRS from park geometries (typically Web Mercator)
            src.crs,  # Target CRS from raster
            always_xy=True
        )
        
        for idx, row in tqdm(parks_gdf.iterrows(), total=len(parks_gdf), desc=f"Calculating {demand_field_prefix} demand metrics"):
            metrics = calculate_demand_metrics_2sfca(
                isochrone_polygon=row["isochrone_polygon"],
                park_geometry=row["geometry"],
                raster=src,
                isochrone_transformer=isochrone_transformer,
                park_transformer=park_transformer,
                decay_type=decay_type
            )
            
            # Add each metric as a column with the prefix
            for metric_name, value in metrics.items():
                col_name = f"{demand_field_prefix}_{metric_name}"
                parks_gdf.at[idx, col_name] = value
    
    return parks_gdf

def demand_analysis(development_polygon: Polygon, population_change: float) -> float:

    # get isochrone polygon for the development polygon
    isochrone_polygon_geojson = fetch_isochrone(development_polygon.centroid.x, development_polygon.centroid.y)
    isochrone_polygon = shape(isochrone_polygon_geojson["geometry"])
    # search radius is the bbox of the isochrone polygon
    # Get radius in degrees, convert to meters
    radius_degrees = isochrone_polygon.centroid.hausdorff_distance(isochrone_polygon.boundary)
    radius_meters = radius_degrees * 111000
    print(f"Searching for parks within {radius_meters} meters of the development polygon at {isochrone_polygon.centroid.x}, {isochrone_polygon.centroid.y}")
    
    parks_gdf = get_parks_gdf_via_arc(
        isochrone_polygon.centroid,  # Already in WGS84
        radius_meters,               # Still need meters
        point_crs="EPSG:4326"        # Tell it the point is in WGS84
    )
    parks_gdf = append_isochrones_to_parks_gdf(parks_gdf)
    parks_gdf = append_demand_metrics(parks_gdf, DEFAULT_GHS_RASTER, "current")

    # Set up the new demand raster
    # clip the population raster to the bounds of all our isochrones
    all_isochrone_bounds = parks_gdf["isochrone_polygon"].total_bounds
    all_isochrone_bounds_polygon = box(*all_isochrone_bounds)
    
    # Get the isochrone CRS (WGS84) - development_polygon is also in this CRS
    isochrone_crs = parks_gdf["isochrone_polygon"].crs
    
    with rasterio.open(DEFAULT_GHS_RASTER) as population_raster:
        # Create transformer: isochrone CRS → raster CRS
        transformer = pyproj.Transformer.from_crs(
            isochrone_crs,  # Source CRS (WGS84)
            population_raster.crs,  # Target CRS from raster
            always_xy=True
        )
        raster_crs = population_raster.crs  # Save for memory raster
        raster_nodata = population_raster.nodata
        
        # Clip population raster to study area (B)
        population_array, _, population_transform = get_raster_clip_under_polygon(
            all_isochrone_bounds_polygon, 
            population_raster, 
            transformer
        )
        
        # Distribute new population in development area (C)
        development_array, development_transform = distribute_population_stats(
            development_polygon, 
            population_change, 
            population_raster, 
            transformer
        )

    # Add the population_clip and the development_clip (with alignment)
    projected_population_array = overlay_rasters(
        population_array, population_transform,
        development_array, development_transform
    )

    # Calculate projected demand using in-memory raster
    with create_memory_raster(projected_population_array, population_transform, raster_crs, raster_nodata) as projected_raster:
        parks_gdf = append_demand_metrics(parks_gdf, projected_raster, "projected")

demand_analysis(development_polygon, 5713)