import rasterio
from rasterio import Affine
from rasterio.io import MemoryFile
import numpy as np
from shapely import Polygon
from shapely.ops import transform
from rasterio.mask import mask
import pyproj
from contextlib import contextmanager

def get_raster_clip_under_polygon(
    polygon: Polygon, 
    raster: rasterio.DatasetReader, 
    transformer=None,
    all_touched: bool = True
) -> tuple[np.ndarray, np.ndarray, Affine]:
    """
    Clip a raster to a polygon and return the clipped array with a valid data mask.
    
    Args:
        polygon: Shapely polygon geometry
        raster: Open rasterio dataset
        transformer: Optional pyproj Transformer to reproject polygon to raster CRS
        all_touched: If True, include all pixels touched by polygon
        
    Returns:
        Tuple of (clipped_array, valid_mask, transform) where:
        - clipped_array: The clipped raster data
        - valid_mask: Boolean mask indicating non-nodata pixels
        - transform: Affine transform for the clipped raster
    """
    # Reproject polygon if transformer provided
    if transformer:
        polygon = transform(transformer.transform, polygon)
    
    out_image, out_transform = mask(
        raster, [polygon], crop=True, all_touched=all_touched
    )
    
    raster_array = out_image[0]
    valid_mask = raster_array != raster.nodata
    
    return raster_array, valid_mask, out_transform

def distribute_population_stats(input_polygon: Polygon, population_change: float, population_raster: rasterio.DatasetReader, transformer: pyproj.Transformer) -> tuple[np.ndarray, rasterio.Affine]:
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


def overlay_rasters(
    base_array: np.ndarray, 
    base_transform: Affine,
    overlay_array: np.ndarray, 
    overlay_transform: Affine
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
    
    # Create result array
    result = base_array.copy().astype(np.float32)
    h, w = overlay_array.shape
    
    # Add overlay at that position (with bounds checking)
    row_end = min(row + h, result.shape[0])
    col_end = min(col + w, result.shape[1])
    row_start = max(row, 0)
    col_start = max(col, 0)
    
    # Calculate overlay slice offsets if row/col are negative
    overlay_row_start = max(0, -row)
    overlay_col_start = max(0, -col)
    overlay_row_end = overlay_row_start + (row_end - row_start)
    overlay_col_end = overlay_col_start + (col_end - col_start)
    
    result[row_start:row_end, col_start:col_end] += overlay_array[
        overlay_row_start:overlay_row_end, 
        overlay_col_start:overlay_col_end
    ]
    
    return result


@contextmanager
def create_memory_raster(
    array: np.ndarray, 
    transform: Affine, 
    crs, 
    nodata=None
):
    """
    Create an in-memory raster from a numpy array that can be read from.
    
    Args:
        array: 2D numpy array of raster data
        transform: Affine transform for the raster
        crs: Coordinate reference system (pyproj CRS or string like "EPSG:4326")
        nodata: Optional nodata value
        
    Yields:
        Open rasterio DatasetReader that can be used with other functions
        
    Example:
        with create_memory_raster(my_array, my_transform, "EPSG:4326") as src:
            # Use src like any rasterio dataset
            data = src.read(1)
    """
    with MemoryFile() as memfile:
        # Write data to memory file
        with memfile.open(
            driver='GTiff',
            height=array.shape[0],
            width=array.shape[1],
            count=1,
            dtype=array.dtype,
            crs=crs,
            transform=transform,
            nodata=nodata
        ) as dst:
            dst.write(array, 1)
        
        # Reopen for reading and yield
        with memfile.open() as src:
            yield src