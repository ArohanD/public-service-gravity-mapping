# main.py
#
# Author: Arohan Dutt
# Date: December 2025
#
# Purpose: Analyze the impact of new residential developments on park demand
#          using a variety of demand metrics. Two-Step Floating Catchment Area
#          (2SFCA) method with distance decay weighting is used to calculate
#          population and acreage demand metrics for each park. See
#          PopulationRaster.py more details on the theory.
#
#
# Procedure Summary:
#          --Take development polygons with expected population changes
#            and a global population raster (GHS-POP).
#          --Fetch driving isochrones from Mapbox API to define
#            catchment areas for developments and nearby parks.
#            Driving is used as this polygon encompasses the other modes
#            of transport. See api.py for more details.
#          --Query OpenStreetMap parks data via ArcGIS feature service
#            for parks that intersect the development catchment areas.
#          --Calculate current demand metrics (population, demand, etc.)
#            using the base population raster.
#          --Distribute new population across rasters clipped to the
#            development polygons.
#          --Recalculate projected demand metrics and calculate delta metrics.
#          --Generate an HTML report with Leaflet maps showing the
#            top impacted parks and their demand changes.
#
#
# Main Steps:
# Step 1: Fetch isochrones for development polygons to define search area.
# Step 2: Query parks that intersect the combined development isochrones.
# Step 3: Fetch isochrones for each park (travel time catchment areas).
# Step 4: Calculate current 2SFCA demand metrics using base population.
# Step 5: Distribute development population onto the raster at development polygon locations.
# Step 6: Calculate projected 2SFCA demand metrics.
# Step 7: Calculate delta metrics (change from current to projected).
# Step 8: Generate HTML report and export shapefile.
#
#
# Inputs:
#   - Development polygons (defined in constants.py as GeoJSON)
#   - Population raster (Optional, default is GHS-POP 100m resolution)
#   - Mapbox API token (for isochrone requests, placed in config.py)
#   - ArcGIS feature service URL (for park polygons)
#
# Outputs:
#   - HTML report with interactive Leaflet maps (park_demand_report.html)
#   - Shapefile with park polygons and demand metrics (placed on an active ArcGIS Pro map)
#
# Software Requirements:
#   - ArcGIS Pro with arcpy
#   - Python packages: geopandas, rasterio, shapely, tqdm, jinja2
#

import os
import arcpy

from utils.html import generate_park_report

# Enable GDAL memory datasets (required for rasterio.mask in newer GDAL versions)
os.environ["GDAL_MEM_ENABLE_OPEN"] = "YES"

from shapely import Polygon
from shapely.geometry import box, shape
from utils.api import fetch_isochrone, get_parks_gdf_via_arc_polygon
from constants import DEVELOPMENT_POLYGONS_GEOJSON
import geopandas as gpd
from tqdm import tqdm
import numpy as np
from utils.arc_utils import gdf_to_featureclass
from utils.utils import log, overlay_rasters, create_memory_raster
from classes.PopulationRaster import PopulationRaster
from shapely.ops import unary_union


DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"

development_polygon_one = Polygon(
    DEVELOPMENT_POLYGONS_GEOJSON["features"][0]["geometry"]["coordinates"][0]
)
development_polygon_two = Polygon(
    DEVELOPMENT_POLYGONS_GEOJSON["features"][1]["geometry"]["coordinates"][0]
)


def append_isochrones_to_parks_gdf(parks_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Fetch isochrones for each park and add as a new column.

    Args:
        parks_gdf: GeoDataFrame with park geometries

    Returns:
        GeoDataFrame with new 'isochrone_polygon' column containing catchment areas
    """
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
    pop_raster: PopulationRaster,
    demand_field_prefix: str,
    decay_type: str = "gaussian",
) -> gpd.GeoDataFrame:
    """Calculate 2SFCA demand metrics for each park's isochrone.

    Args:
        parks_gdf: GeoDataFrame with park geometries and isochrone_polygon column
        pop_raster: PopulationRaster instance to calculate demand from
        demand_field_prefix: Prefix for output columns (e.g., "current" -> "current_pop", etc.)
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
    isochrone_crs = parks_gdf["isochrone_polygon"].crs
    park_crs = parks_gdf.crs

    for idx, row in tqdm(
        parks_gdf.iterrows(),
        total=len(parks_gdf),
        desc=f"Calculating {demand_field_prefix} demand metrics",
    ):
        metrics = pop_raster.calculate_demand_metrics(
            isochrone_polygon=row["isochrone_polygon"],
            park_geometry=row["geometry"],
            isochrone_crs=isochrone_crs,
            park_crs=park_crs,
            decay_type=decay_type,
        )

        # Add each metric as a column with the prefix
        for metric_name, value in metrics.items():
            col_name = f"{demand_field_prefix}_{metric_name}"
            parks_gdf.at[idx, col_name] = value

    return parks_gdf


def append_delta_metrics(parks_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Calculate change metrics between current and projected demand.

    Requires parks_gdf to have current_* and projected_* columns from append_demand_metrics.

    Output columns added:
        pop_change, m2_per_person_change, acres_per_1000_change: Raw differences
        pop_weighted_change, m2_per_person_weighted_change, acres_per_1000_weighted_change: Weighted differences
        pop_pct_change, acres_per_1000_pct_change: Percent changes (raw)
        pop_weighted_pct_change, acres_per_1000_weighted_pct_change: Percent changes (weighted)
    """
    # Calculate change from current to projected (positive = increase)
    # Raw change stats
    parks_gdf["pop_change"] = parks_gdf["projected_pop"] - parks_gdf["current_pop"]
    parks_gdf["m2_per_person_change"] = (
        parks_gdf["projected_m2_per_person"] - parks_gdf["current_m2_per_person"]
    )
    parks_gdf["acres_per_1000_change"] = (
        parks_gdf["projected_acres_per_1000"] - parks_gdf["current_acres_per_1000"]
    )

    # Weighted change stats (distance-decay weighted)
    parks_gdf["pop_weighted_change"] = (
        parks_gdf["projected_pop_weighted"] - parks_gdf["current_pop_weighted"]
    )
    parks_gdf["m2_per_person_weighted_change"] = (
        parks_gdf["projected_m2_per_person_weighted"]
        - parks_gdf["current_m2_per_person_weighted"]
    )
    parks_gdf["acres_per_1000_weighted_change"] = (
        parks_gdf["projected_acres_per_1000_weighted"]
        - parks_gdf["current_acres_per_1000_weighted"]
    )

    # Percent change (useful for comparing across different sized parks)
    parks_gdf["pop_pct_change"] = (
        parks_gdf["pop_change"] / parks_gdf["current_pop"]
    ) * 100
    parks_gdf["acres_per_1000_pct_change"] = (
        parks_gdf["acres_per_1000_change"] / parks_gdf["current_acres_per_1000"]
    ) * 100

    # Weighted percent change
    parks_gdf["pop_weighted_pct_change"] = (
        parks_gdf["pop_weighted_change"] / parks_gdf["current_pop_weighted"]
    ) * 100
    parks_gdf["acres_per_1000_weighted_pct_change"] = (
        parks_gdf["acres_per_1000_weighted_change"]
        / parks_gdf["current_acres_per_1000_weighted"]
    ) * 100

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
        log("Warning: developments_gdf has no CRS, assuming WGS84")
        developments_gdf = developments_gdf.set_crs("EPSG:4326")
    elif developments_gdf.crs != "EPSG:4326":
        developments_gdf = developments_gdf.to_crs("EPSG:4326")

    # ------------------
    # Step 1: Fetch isochrones for development polygons to define search area.
    log(f"Processing {len(developments_gdf)} development polygons...")

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

    # ------------------
    # Step 2: Query parks that intersect the combined development isochrones.
    combined_isochrone = unary_union(development_isochrones)
    log("Searching for parks that intersect the combined development isochrone area")

    parks_gdf = get_parks_gdf_via_arc_polygon(
        combined_isochrone, polygon_crs="EPSG:4326"
    )

    if len(parks_gdf) == 0:
        log("No parks found in the search area")
        return parks_gdf

    # ------------------
    # Step 3: Fetch isochrones for each park (travel time catchment areas).
    parks_gdf = append_isochrones_to_parks_gdf(parks_gdf)

    # Open the population raster once for the entire analysis
    with PopulationRaster(DEFAULT_GHS_RASTER) as pop_raster:
        # ------------------
        # Step 4: Calculate current 2SFCA demand metrics using base population.
        parks_gdf = append_demand_metrics(parks_gdf, pop_raster, "current")

        # ------------------
        # Step 5: Distribute development population onto the raster.
        study_area = unary_union(
            [parks_gdf["isochrone_polygon"].union_all(), developments_gdf.union_all()]
        )
        all_isochrone_bounds_polygon = box(*study_area.bounds)
        isochrone_crs = parks_gdf["isochrone_polygon"].crs

        # Clip population raster to study area (base layer)
        population_array, _, population_transform = pop_raster.clip_to_polygon(
            all_isochrone_bounds_polygon, polygon_crs=isochrone_crs
        )

        # Start with the base population
        projected_population_array = population_array.copy().astype(np.float32)

        # Overlay each development polygon's population onto the base
        for _, row in tqdm(
            developments_gdf.iterrows(),
            total=len(developments_gdf),
            desc="Distributing development population",
        ):
            development_polygon = row.geometry
            population_change = row["population_change"]

            if population_change == 0:
                continue

            # Distribute population for this development
            development_array, development_transform = pop_raster.distribute_population(
                development_polygon, population_change, polygon_crs=isochrone_crs
            )

            # Overlay onto the projected population
            projected_population_array = overlay_rasters(
                projected_population_array,
                population_transform,
                development_array,
                development_transform,
            )

        # ------------------
        # Step 6: Calculate projected 2SFCA demand metrics.
        with create_memory_raster(
            projected_population_array,
            population_transform,
            pop_raster.crs,
            pop_raster.nodata,
        ) as projected_raster:
            parks_gdf = append_demand_metrics(parks_gdf, projected_raster, "projected")

    # ------------------
    # Step 7: Calculate delta metrics (change from current to projected).
    parks_gdf = append_delta_metrics(parks_gdf)

    return parks_gdf


# Create a GeoDataFrame with development polygons and their population changes
developments_gdf = gpd.GeoDataFrame(
    {
        "geometry": [development_polygon_one, development_polygon_two],
        "population_change": [500.0, 1000.0],
    },
    crs="EPSG:4326",
)

# Run the analysis
arcpy.AddMessage("Starting demand analysis...")
parks_gdf = demand_analysis(developments_gdf)

# ------------------
# Step 8: Generate HTML report and export shapefile.
html_report_path = generate_park_report(parks_gdf)
os.startfile(html_report_path)

arcpy.AddMessage(f"Analysis complete. Adding {len(parks_gdf)} parks to map...")

# Save to a shapefile instead of in-memory (more reliable)
output_shp = r"C:\gispy\scratch\parks_demand_analysis.shp"

fc_path = gdf_to_featureclass(parks_gdf, output_path=output_shp)
arcpy.AddMessage(f"Saved to: {fc_path}")

# Add to map (only works when running inside ArcGIS Pro)
try:
    aprx = arcpy.mp.ArcGISProject("CURRENT")
    m = aprx.activeMap
    if m:
        arcpy.AddMessage(f"Active map: {m.name}")
        m.addDataFromPath(fc_path)
        arcpy.AddMessage("Layer added!")
    else:
        arcpy.AddWarning("No active map found")
except Exception:
    # "CURRENT" error means we're running from command line, not ArcGIS Pro
    log(f"Note: Not running in ArcGIS Pro - shapefile saved to: {fc_path}")
    log("You can manually add this shapefile to your map.")
