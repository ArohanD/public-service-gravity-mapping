# inputs: development polygon, population raster
#


from shapely import Polygon, Point
from get_parks import get_parks_gdf


DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"


def demand_analysis(development_polygon: Polygon, population_raster: str) -> float:
    pass


# Test parks gdf fetch
parks_gdf = get_parks_gdf(Point(-8786562.4714, 4300726.410700001), 50)
import pdb; pdb.set_trace()