import arcpy
import rasterio
import os
import numpy as np
from shapely import Polygon, wkt
import pyproj

from arc_utils import update_map
from utils import get_raster_clip_under_polygon

# Enable GDAL memory datasets (required for rasterio.mask in newer GDAL versions)
os.environ["GDAL_MEM_ENABLE_OPEN"] = "YES"

# Script tool params
isochrones_layer = arcpy.GetParameter(0)
ghs_raster = arcpy.GetParameterAsText(1)
out_fc = arcpy.GetParameterAsText(2)

arcpy.env.overwriteOutput = True

# Default params for standalone execution
DEFAULT_ISCHRONE_LAYER = r"..\Data\samples\park_isochrones_1000.shp"  # TODO: Change this to the actual isochrone layer
DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"
DEFAULT_OUT_FC = r"..\gispy\scratch\park_isochrones_with_demand.shp"

DEMAND_FIELD_NAME = "DEMAND"


def add_demand_field(fc, field_name):
    """Add the demand field to store population sums."""
    existing_fields = [f.name for f in arcpy.ListFields(fc)]
    if field_name not in existing_fields:
        arcpy.management.AddField(fc, field_name, "DOUBLE")
    arcpy.AddMessage(f"Added field: {field_name}")


def calculate_demand_for_polygon(polygon: Polygon, raster: rasterio.DatasetReader, transformer):
    """
    Calculate demand for a polygon by summing the population within the polygon from the raster.
    
    Args:
        polygon: Shapely polygon geometry
        raster: Open rasterio dataset
        transformer: pyproj Transformer to reproject polygon to raster CRS
    """
    try:
        raster_array, valid_mask, _ = get_raster_clip_under_polygon(polygon, raster, transformer, all_touched=True)
        demand = float(np.sum(raster_array[valid_mask]))
        return demand

    except Exception as e:
        arcpy.AddError(f"Error calculating demand for polygon: {e}")
        return None


def calculate_demand(isochrones, raster_path, out_fc):
    """Calculate demand at each isochrone polygon by summing the population within the polygon from the raster."""
    arcpy.AddMessage(f"Calculating demand under {isochrones} using {raster_path}")

    # copy isochrones to output feature class
    if not out_fc:
        out_fc = DEFAULT_OUT_FC
    arcpy.management.CopyFeatures(isochrones, out_fc)

    add_demand_field(out_fc, DEMAND_FIELD_NAME)

    total_isochrones = arcpy.management.GetCount(out_fc)[0]
    arcpy.AddMessage(f"Processing {total_isochrones} isochrones...")

    # Get source CRS from isochrones
    desc = arcpy.Describe(out_fc)
    source_crs = desc.spatialReference
    arcpy.AddMessage(f"Isochrones CRS: {source_crs.name} (WKID: {source_crs.factoryCode})")

    with rasterio.open(raster_path) as src:
        raster_crs = src.crs
        arcpy.AddMessage(f"Raster CRS: {raster_crs}")
        arcpy.AddMessage(f"Raster bounds: {src.bounds}")
        
        # Create transformer from isochrones CRS to raster CRS
        transformer = pyproj.Transformer.from_crs(
            f"EPSG:{source_crs.factoryCode}",  # Source CRS from isochrones
            raster_crs.to_string(),  # Target CRS from raster
            always_xy=True
        )
        
        # Keep track for progress updates
        processed_count = 0

        with arcpy.da.UpdateCursor(out_fc, ["SHAPE@", DEMAND_FIELD_NAME]) as cursor:
            for row in cursor:
                polygon = row[0]
                shapely_polygon = wkt.loads(polygon.WKT)  # Convert to shapely polygon
                demand = calculate_demand_for_polygon(shapely_polygon, src, transformer)

                row[1] = demand
                cursor.updateRow(row)
                processed_count += 1
                arcpy.AddMessage(
                    f"Processed {processed_count} of {total_isochrones} isochrones"
                )

    arcpy.AddMessage("Demand calculation complete.")

    update_map(out_fc)
    return out_fc


if not isochrones_layer:
    isochrones = DEFAULT_ISCHRONE_LAYER
else:
    isochrones = isochrones_layer

if not ghs_raster:
    raster_path = DEFAULT_GHS_RASTER
else:
    raster_path = ghs_raster

if not out_fc:
    out_fc = DEFAULT_OUT_FC

calculate_demand(isochrones, raster_path, out_fc)
