"""Population raster utilities for demand analysis."""

from rasterio import Affine
from rasterio.io import MemoryFile
import numpy as np
from contextlib import contextmanager
from classes.PopulationRaster import PopulationRaster

def overlay_rasters(
    base_array: np.ndarray,
    base_transform: Affine,
    overlay_array: np.ndarray,
    overlay_transform: Affine,
) -> np.ndarray:
    """
    Add a smaller overlay raster onto a larger base raster at the correct position.
    
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


def create_population_raster_from_array(
    array: np.ndarray,
    transform: Affine,
    crs,
    nodata=None
) -> "PopulationRaster":
    """
    Create a PopulationRaster from a numpy array (in-memory).
    
    Args:
        array: 2D numpy array of population data
        transform: Affine transform for the raster
        crs: Coordinate reference system
        nodata: Optional nodata value
        
    Returns:
        PopulationRaster wrapping the in-memory data
        
    Note: The returned PopulationRaster must be used within a MemoryFile context.
    Use create_memory_raster() for a simpler context manager approach.
    """
    memfile = MemoryFile()
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
    
    src = memfile.open()
    raster = PopulationRaster(src)
    raster._memfile = memfile  # Keep reference to prevent garbage collection
    return raster


# Keep legacy function for backward compatibility
def create_memory_raster(array: np.ndarray, transform: Affine, crs, nodata=None):
    """
    Context manager to create an in-memory raster from a numpy array.
    
    Yields a PopulationRaster that can be used for demand calculations.
    
    Example:
        with create_memory_raster(my_array, my_transform, "EPSG:4326") as pop_raster:
            metrics = pop_raster.calculate_demand_metrics(...)
    """
    
    @contextmanager
    def _create():
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
    
    return _create()
