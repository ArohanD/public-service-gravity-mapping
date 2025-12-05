import arcpy
import json
import urllib.request
import urllib.parse
from shapely import Point
import geopandas as gpd


from config import MAPBOX_TOKEN

arcpy.env.overwriteOutput = True


def fetch_isochrone(lon, lat, mode="driving-traffic", token=MAPBOX_TOKEN) -> dict:
    """
    Fetch an isochrone from the Mapbox API for a given longitude and latitude.
    Returns a JSON object containing the isochrone polygons, or None on error.
    """
    mode = "walking"  # reducing range for now to reduce API calls TODO: add back in driving-traffic
    base_url = f"https://api.mapbox.com/isochrone/v1/mapbox/{mode}"
    coords = f"{lon},{lat}"

    params = {
        "contours_minutes": "10",
        "polygons": "true",
        "denoise": "1",
        "access_token": token,
    }

    url = f"{base_url}/{coords}?{urllib.parse.urlencode(params)}"

    try:
        with urllib.request.urlopen(url) as response:
            text = response.read().decode("utf-8")
        parsed_response = json.loads(text)
        feature = parsed_response["features"][0]
        return feature
    except Exception as e:
        arcpy.AddError(f"Failed to fetch isochrone for ({lon}, {lat}): {e}")
        return None


# Get parks within a radius of a point using ArcGIS feature services.

# Default parameters
DEFAULT_WEB_LAYER = "https://services6.arcgis.com/Do88DoK2xjTUCXd1/arcgis/rest/services/OSM_NA_Leisure/FeatureServer/0"
DEFAULT_DEFINITION_QUERY = "leisure = 'park'"


def get_parks_gdf_via_arc(
    point: Point,
    radius: float,
    parks: str = DEFAULT_WEB_LAYER,
    point_crs: str = "EPSG:4326",
) -> gpd.GeoDataFrame:
    """Get parks within a radius of a point using arcpy. Returns a GeoDataFrame.

    Args:
        point: Shapely Point geometry
        radius: Buffer radius in meters
        parks: URL to parks feature service
        point_crs: CRS of the input point (default: EPSG:4326 for WGS84)
    """
    from shapely import wkt

    # Clean up existing layers
    if arcpy.Exists("parks_lyr"):
        arcpy.management.Delete("parks_lyr")
    if arcpy.Exists("in_memory\\click_buffer"):
        arcpy.management.Delete("in_memory\\click_buffer")

    # Get EPSG code from point_crs
    epsg_code = int(point_crs.split(":")[1])
    sr = arcpy.SpatialReference(epsg_code)

    # Create point geometry and buffer
    pt_geom = arcpy.PointGeometry(arcpy.Point(point.x, point.y), sr)
    buffer_fc = "in_memory\\click_buffer"
    arcpy.analysis.Buffer(pt_geom, buffer_fc, f"{radius} Meters")

    # Create parks layer with definition query
    arcpy.management.MakeFeatureLayer(parks, "parks_lyr", DEFAULT_DEFINITION_QUERY)

    # Get the parks layer's CRS (geometries will be in this CRS)
    parks_sr = arcpy.Describe("parks_lyr").spatialReference
    parks_epsg = parks_sr.factoryCode

    # Select parks that intersect buffer
    arcpy.management.SelectLayerByLocation("parks_lyr", "INTERSECT", buffer_fc)

    # Get count
    count = int(arcpy.management.GetCount("parks_lyr")[0])
    print(f"Found {count} park polygons.")

    if count == 0:
        return gpd.GeoDataFrame()

    # Convert selected features to GeoDataFrame
    geometries = []
    attributes = []

    # Get field names (excluding geometry and OID fields)
    fields = [
        f.name
        for f in arcpy.ListFields("parks_lyr")
        if f.type not in ("Geometry", "OID")
    ]

    with arcpy.da.SearchCursor("parks_lyr", ["SHAPE@WKT"] + fields) as cursor:
        for row in cursor:
            geom = wkt.loads(row[0])
            geometries.append(geom)
            attributes.append(dict(zip(fields, row[1:])))

    # Create GeoDataFrame with the parks layer's CRS (not the input point's CRS)
    gdf = gpd.GeoDataFrame(attributes, geometry=geometries, crs=f"EPSG:{parks_epsg}")

    return gdf
