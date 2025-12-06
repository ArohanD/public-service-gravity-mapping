# utils/utils.py
#
# Author: Arohan Dutt
# Date: December 2025
#
# Purpose: Various utility functions for the project.
#

from rasterio import Affine
from rasterio.io import MemoryFile
import numpy as np
from contextlib import contextmanager
from classes.PopulationRaster import PopulationRaster
import arcpy


def log(*args, **kwargs):
    # if arcpy, use arcpy.AddMessage
    try:
        arcpy.mp.ArcGISProject("CURRENT")
        arcpy.AddMessage(*args, **kwargs)
    except OSError:
        print(*args, **kwargs)


def overlay_rasters(
    base_array: np.ndarray,
    base_transform: Affine,
    overlay_array: np.ndarray,
    overlay_transform: Affine,
) -> np.ndarray:
    """
    Add a smaller overlay raster onto a larger base raster at the correct position.
    Adapted from:
    https://gis.stackexchange.com/questions/250077/overlapping-rasters-as-numpy-arrays

    Args:
        base_array: The base raster array (larger)
        base_transform: Affine transform for the base raster
        overlay_array: The overlay raster array (smaller)
        overlay_transform: Affine transform for the overlay raster

    Returns:
        Combined array with overlay added to base at the correct position
    """
    # Get overlay's top-left corner in world coordinates
    overlay_x = overlay_transform.c
    overlay_y = overlay_transform.f

    # Convert to base's pixel coordinates using inverse transform
    col, row = ~base_transform * (overlay_x, overlay_y)
    col, row = int(round(col)), int(round(row))

    result = base_array.copy().astype(np.float32)
    h, w = overlay_array.shape

    # Calculate slice bounds with bounds checking
    row_end = min(row + h, result.shape[0])
    col_end = min(col + w, result.shape[1])
    row_start = max(row, 0)
    col_start = max(col, 0)

    overlay_row_start = max(0, -row)
    overlay_col_start = max(0, -col)
    overlay_row_end = overlay_row_start + (row_end - row_start)
    overlay_col_end = overlay_col_start + (col_end - col_start)

    result[row_start:row_end, col_start:col_end] += overlay_array[
        overlay_row_start:overlay_row_end, overlay_col_start:overlay_col_end
    ]

    return result


@contextmanager
def create_memory_raster(array: np.ndarray, transform: Affine, crs, nodata=None):
    """
    Context manager to create an in-memory PopulationRaster from a numpy array.
    Adapted from: https://gist.github.com/lpinner/13244b5c589cda4fbdfa89b30a44005b

    Args:
        array: 2D numpy array of population data
        transform: Affine transform for the raster
        crs: Coordinate reference system
        nodata: Optional nodata value

    Yields:
        PopulationRaster instance

    Example:
        with create_memory_raster(projected_array, transform, crs) as pop_raster:
            parks_gdf = append_demand_metrics(parks_gdf, pop_raster, "projected")
    """
    with MemoryFile() as memfile:
        with memfile.open(
            driver="GTiff",
            height=array.shape[0],
            width=array.shape[1],
            count=1,
            dtype=array.dtype,
            crs=crs,
            transform=transform,
            nodata=nodata,
        ) as dst:
            dst.write(array, 1)

        with memfile.open() as src:
            yield PopulationRaster(src)
