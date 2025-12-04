import rasterio
import numpy as np
from shapely import Polygon
from shapely.ops import transform
from rasterio.mask import mask

def get_raster_clip_under_polygon(
    polygon: Polygon, 
    raster: rasterio.DatasetReader, 
    transformer=None,
    all_touched: bool = True
) -> tuple[np.ndarray, np.ndarray]:
    """
    Clip a raster to a polygon and return the clipped array with a valid data mask.
    
    Args:
        polygon: Shapely polygon geometry
        raster: Open rasterio dataset
        transformer: Optional pyproj Transformer to reproject polygon to raster CRS
        all_touched: If True, include all pixels touched by polygon
        
    Returns:
        Tuple of (clipped_array, valid_mask) where valid_mask indicates non-nodata pixels
    """
    # Reproject polygon if transformer provided
    if transformer:
        polygon = transform(transformer.transform, polygon)
    
    out_image, _out_transform = mask(
        raster, [polygon], crop=True, all_touched=all_touched
    )
    
    raster_array = out_image[0]
    valid_mask = raster_array != raster.nodata
    
    return raster_array, valid_mask