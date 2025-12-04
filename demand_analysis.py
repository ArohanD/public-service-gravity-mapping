# inputs: development polygon, population raster
#


from shapely import Polygon, Point
from shapely.geometry import shape
from get_isochrones import fetch_isochrone
from get_parks import get_parks_gdf_via_arc
from constants import DEVELOPMENT_POLYGON_GEOJSON

DEFAULT_GHS_RASTER = r"..\Data\GHS\GHS_POP_E2025_GLOBE_R2023A_54009_100_V1_0.tif"

development_polygon = Polygon(DEVELOPMENT_POLYGON_GEOJSON["features"][0]["geometry"]["coordinates"][0])


def demand_analysis(development_polygon: Polygon, population_change: float) -> float:

    # get isochrone polygon for the development polygon
    isochrone_polygon_geojson = fetch_isochrone(development_polygon.centroid.x, development_polygon.centroid.y)
    isochrone_polygon = shape(isochrone_polygon_geojson["geometry"])
    # search radius is the bbox of the isochrone polygon
    # Get radius in degrees, convert to meters
    radius_degrees = isochrone_polygon.centroid.hausdorff_distance(isochrone_polygon.boundary)
    radius_meters = radius_degrees * 111000
    print(f"Searching for parks within {radius_meters} meters of the development polygon at {isochrone_polygon.centroid.x}, {isochrone_polygon.centroid.y}")
    
    parks_gdf = get_parks_gdf_via_arc(
        isochrone_polygon.centroid,  # Already in WGS84
        radius_meters,               # Still need meters
        point_crs="EPSG:4326"        # Tell it the point is in WGS84
    )
    import pdb; pdb.set_trace()


# Test parks gdf fetch
# parks_gdf = get_parks_gdf(Point(-8786562.4714, 4300726.410700001), 50)

demand_analysis(development_polygon, 5713)