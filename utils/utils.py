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

def apply_distance_decay(
    pop_array: np.ndarray,
    transform,
    target_x: float,
    target_y: float,
    decay_type: str = "gaussian",
    sigma_km: float = 0.5
) -> np.ndarray:
    """
    Weight population cells by distance to a target point.
    
    Args:
        pop_array: Population values per cell
        transform: Raster affine transform (gives cell coordinates)
        target_x/y: Target location in raster CRS (e.g., park centroid)
        decay_type: "none", "inverse", "inverse_square", or "gaussian"
        sigma_km: For gaussian decay, controls spread (higher = slower decay)
    
    Returns:
        Weighted population array
    """
    if decay_type == "none":
        return pop_array
    
    rows, cols = pop_array.shape
    
    # Create coordinate grids for all cells
    col_indices, row_indices = np.meshgrid(np.arange(cols), np.arange(rows))
    
    # Convert pixel indices to world coordinates (cell centers)
    cell_x = transform.c + (col_indices + 0.5) * transform.a
    cell_y = transform.f + (row_indices + 0.5) * transform.e
    
    # Calculate distance from each cell to target (in raster CRS units, typically meters)
    distances = np.sqrt((cell_x - target_x)**2 + (cell_y - target_y)**2)
    distances_km = distances / 1000
    
    # Apply decay function
    if decay_type == "inverse":
        weights = 1 / (1 + distances_km)
    elif decay_type == "inverse_square":
        weights = 1 / (1 + distances_km**2)
    elif decay_type == "gaussian":
        weights = np.exp(-distances_km**2 / (2 * sigma_km**2))
    else:
        weights = np.ones_like(distances)
    
    return pop_array * weights


def calculate_demand_metrics_2sfca(
    isochrone_polygon: Polygon,
    park_geometry: Polygon,
    raster: rasterio.DatasetReader,
    isochrone_transformer: pyproj.Transformer,
    park_transformer: pyproj.Transformer,
    decay_type: str = "gaussian"
) -> dict:
    """
    Calculate 2SFCA demand metrics for a park.
    
    Args:
        isochrone_polygon: The catchment area (isochrone) for the park
        park_geometry: The park polygon (for area calculation)
        raster: Open rasterio dataset
        isochrone_transformer: Transformer from isochrone CRS to raster CRS
        park_transformer: Transformer from park geometry CRS to raster CRS
        decay_type: Distance decay function type
    
    Returns:
        Dictionary with multiple demand metrics
    """
    from shapely.ops import transform as shapely_transform
    
    try:
        # Get population within isochrone
        raster_array, valid_mask, out_transform = get_raster_clip_under_polygon(
            isochrone_polygon, raster, isochrone_transformer, all_touched=True
        )
        
        # Transform park centroid to raster CRS for distance calculations
        # (park geometry is in different CRS than isochrone)
        park_centroid = park_geometry.centroid
        park_centroid_raster = shapely_transform(park_transformer.transform, park_centroid)
        
        # Calculate weighted population (with distance decay)
        weighted_pop = apply_distance_decay(
            raster_array,
            out_transform,
            park_centroid_raster.x,
            park_centroid_raster.y,
            decay_type=decay_type
        )
        
        # Raw population sum (no decay)
        pop_raw = float(np.sum(raster_array[valid_mask]))
        
        # Weighted population sum (with decay)
        pop_weighted = float(np.sum(weighted_pop[valid_mask]))
        
        # Calculate park area in m² (transform to raster CRS which is typically in meters)
        park_raster_crs = shapely_transform(park_transformer.transform, park_geometry)
        park_area_m2 = park_raster_crs.area
        
        # 2SFCA metrics
        # Rj = Supply / Demand
        m2_per_person = park_area_m2 / pop_raw if pop_raw > 0 else None
        m2_per_person_weighted = park_area_m2 / pop_weighted if pop_weighted > 0 else None
        
        # More readable units
        acres_per_1000 = (park_area_m2 / 4047) / (pop_raw / 1000) if pop_raw > 0 else None
        acres_per_1000_weighted = (park_area_m2 / 4047) / (pop_weighted / 1000) if pop_weighted > 0 else None
        
        return {
            "pop": pop_raw,
            "pop_weighted": pop_weighted,
            "park_area_m2": park_area_m2,
            "m2_per_person": m2_per_person,
            "m2_per_person_weighted": m2_per_person_weighted,
            "acres_per_1000": acres_per_1000,
            "acres_per_1000_weighted": acres_per_1000_weighted,
        }

    except Exception as e:
        print(f"Error calculating demand metrics: {e}")
        return {
            "pop": None,
            "pop_weighted": None,
            "park_area_m2": None,
            "m2_per_person": None,
            "m2_per_person_weighted": None,
            "acres_per_1000": None,
            "acres_per_1000_weighted": None,
        }