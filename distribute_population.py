# takes an input polygon with a population
# clips the polygon over the population raster
# distributes the population to the cells that are inside the polygon

import arcpy
import os
import rasterio
import pyproj
import sys
from utils.arc_utils import load_first_polygon_from_shapefile
from utils.utils import distribute_population_stats

# Enable GDAL memory datasets (required for rasterio.mask in newer GDAL versions)
os.environ['GDAL_MEM_ENABLE_OPEN'] = 'YES'
arcpy.env.overwriteOutput = True

DEFAULT_INPUT_POLYGON = r"..\Data\samples\developmentPolygon.shp"  # TODO: Change this to the actual single polygon layer
DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"
DEFAULT_OUT_RASTER = r"c:\gispy\scratch\distributed_population_raster.tif"

input_polygon = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT_POLYGON
population_change = float(sys.argv[2]) if len(sys.argv) > 2 else 5713
population_raster = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_GHS_RASTER

if __name__ == "__main__":
    desc = arcpy.Describe(input_polygon)
    polygon_crs = desc.spatialReference.factoryCode

    with rasterio.open(population_raster) as src:
        # load shapfile into a polygon
        input_polygon = load_first_polygon_from_shapefile(input_polygon)
        raster_crs = src.crs
        transformer = pyproj.Transformer.from_crs(
            f"EPSG:{polygon_crs}",  # Source CRS from input polygon
            raster_crs.to_string(),  # Target CRS from raster
            always_xy=True
        )
        population_array, out_transform = distribute_population_stats(input_polygon, population_change, src, transformer)

        # save the output raster with correct dimensions for the clipped extent
        profile = src.profile.copy()
        profile.update(
            height=population_array.shape[0],
            width=population_array.shape[1],
            transform=out_transform
        )
        
        with rasterio.open(DEFAULT_OUT_RASTER, 'w', **profile) as dst:
            dst.write(population_array, 1)