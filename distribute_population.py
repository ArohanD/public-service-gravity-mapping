# takes an input polygon with a population
# clips the polygon over the population raster
# distributes the population to the cells that are inside the polygon

import arcpy
import os
import rasterio
import numpy as np
import pyproj
import sys
from shapely import Polygon
from arc_utils import load_first_polygon_from_shapefile
from utils.utils import get_raster_clip_under_polygon

# Enable GDAL memory datasets (required for rasterio.mask in newer GDAL versions)
os.environ['GDAL_MEM_ENABLE_OPEN'] = 'YES'
arcpy.env.overwriteOutput = True

DEFAULT_INPUT_POLYGON = r"..\Data\samples\developmentPolygon.shp"  # TODO: Change this to the actual single polygon layer
DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"
DEFAULT_OUT_RASTER = r"c:\gispy\scratch\distributed_population_raster.tif"

def distribute_population(input_polygon: Polygon, population_change: float, population_raster: rasterio.DatasetReader, transformer: pyproj.Transformer) -> tuple[np.ndarray, rasterio.Affine]:
    """
    Distribute the population to the cells that are inside the polygon (randomly).
    
    Returns:
        Tuple of (output_array, transform) for writing to a raster file.
    """

    # clip the population raster to the input polygon
    raster_array, valid_mask, out_transform = get_raster_clip_under_polygon(input_polygon, population_raster, transformer)
    # distribute the population to the cells that are inside the polygon

    # Note: The NoData value should only appear in uninhabited areas,
    # cells with no population have a value of 0 and are still candidates
    # for distribution.
    valid_cells = np.where(valid_mask)
    num_valid_cells = len(valid_cells[0])

    if num_valid_cells == 0:
        raise ValueError("No valid cells found in the polygon")

    # Give each current cell a random weight that in adds up to the population change
    random_weights = np.random.rand(num_valid_cells)
    random_weights /= np.sum(random_weights) # normalize to 1
    output = np.zeros(raster_array.shape, dtype=np.float32)
    output[valid_cells] = random_weights * population_change

    # return the output raster and its transform
    return output, out_transform


input_polygon = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_INPUT_POLYGON
population_change = float(sys.argv[2]) if len(sys.argv) > 2 else 5713
population_raster = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_GHS_RASTER

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
    population_array, out_transform = distribute_population(input_polygon, population_change, src, transformer)

    # save the output raster with correct dimensions for the clipped extent
    profile = src.profile.copy()
    profile.update(
        height=population_array.shape[0],
        width=population_array.shape[1],
        transform=out_transform
    )
    
    with rasterio.open(DEFAULT_OUT_RASTER, 'w', **profile) as dst:
        dst.write(population_array, 1)
