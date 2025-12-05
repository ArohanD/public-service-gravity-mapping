# inputs: development polygon, population raster
#


import pyproj
from shapely import Polygon
from shapely.geometry import box, shape
from calculate_demand_at_park_isochrones import calculate_demand_for_polygon
from get_isochrones import fetch_isochrone
from get_parks import get_parks_gdf_via_arc
from constants import DEVELOPMENT_POLYGON_GEOJSON
import geopandas as gpd
from tqdm import tqdm
import rasterio
import numpy as np

from utils.utils import distribute_population_stats, get_raster_clip_under_polygon, overlay_rasters, create_memory_raster

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


def calculate_demand_metrics(polygon: Polygon, raster: rasterio.DatasetReader, transformer: pyproj.Transformer) -> float:
    """Calculate demand for a polygon by summing the population within the polygon from the raster."""
    current_demand = calculate_demand_for_polygon(polygon, raster, transformer)
    
    return current_demand, None


def append_demand_metrics(parks_gdf: gpd.GeoDataFrame, population_raster: str, demand_field_name: str) -> gpd.GeoDataFrame:
    with rasterio.open(population_raster) as src:
        transformer = pyproj.Transformer.from_crs(
            parks_gdf["isochrone_polygon"].crs,  # Source CRS from input polygon
            src.crs,  # Target CRS from raster
            always_xy=True
        )
        for _, row in tqdm(parks_gdf.iterrows(), total=len(parks_gdf), desc="Calculating demand at park isochrones"):
            current_demand, projected_demand = calculate_demand_metrics(row["isochrone_polygon"], src, transformer)
            parks_gdf.at[row.name, demand_field_name] = current_demand
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
    parks_gdf = append_demand_metrics(parks_gdf, DEFAULT_GHS_RASTER, "current_demand")

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
        for _, row in tqdm(parks_gdf.iterrows(), total=len(parks_gdf), desc="Calculating projected demand"):
            demand, _ = calculate_demand_metrics(row["isochrone_polygon"], projected_raster, transformer)
            parks_gdf.at[row.name, "projected_demand"] = demand

demand_analysis(development_polygon, 5713)