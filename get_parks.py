# Input a point and a radius (in meters) and then output polygons
# representing all parks that touch the radius.

import arcpy
import sys
from urllib.parse import quote
from get_point import get_point
from shapely import Point
import geopandas as gpd

arcpy.env.overwriteOutput = True

# Default parameters
DEFAULT_WEB_LAYER = "https://services6.arcgis.com/Do88DoK2xjTUCXd1/arcgis/rest/services/OSM_NA_Leisure/FeatureServer/0"
DEFAULT_DEFINITION_QUERY = "leisure = 'park'"
DEFAULT_RADIUS = 50.0
DEFAULT_COORDS = (-8784459.1687, 4294077.940800004)
DEFAULT_OUT_FC = "c:\\gispy\\scratch\\parks_in_radius.shp"


def parse_cli_args():
    """Parse command line arguments for standalone mode."""
    if len(sys.argv) > 1:
        coords = sys.argv[1].split(',')
        x, y = float(coords[0]), float(coords[1])
    else:
        x, y = DEFAULT_COORDS
    
    radius = float(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_RADIUS
    parks = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_WEB_LAYER
    out_fc = sys.argv[4] if len(sys.argv) > 4 else DEFAULT_OUT_FC
    
    return (x, y), radius, parks, out_fc


def get_parks(point, radius, parks=DEFAULT_WEB_LAYER, out_fc=DEFAULT_OUT_FC):
    """Get parks within a radius of a point. Returns a feature class of parks that intersect the buffer."""


    running_in_arc = not isinstance(point, tuple)
    
    # Clean up existing layers
    if arcpy.Exists("parks_lyr"):
        arcpy.management.Delete("parks_lyr")
    if arcpy.Exists("in_memory\\click_buffer"):
        arcpy.management.Delete("in_memory\\click_buffer")
    
    # Extract coordinates from point (handle both feature and tuple)
    x, y = get_point(point) if running_in_arc else point
    arcpy.AddMessage(f"Getting parks within a radius of {radius} meters from point: {x}, {y}")
    
    # Get spatial reference from the INPUT point (not parks layer) to avoid CRS mismatch
    if running_in_arc:
        sr = arcpy.Describe(point).spatialReference
    else:
        sr = arcpy.Describe(parks).spatialReference
    
    # Create point geometry and buffer in the point's CRS
    pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)
    buffer_fc = "in_memory\\click_buffer"
    arcpy.analysis.Buffer(pt_geom, buffer_fc, f"{radius} Meters")

    # Create parks layer with definition query if using URL
    definition_query = DEFAULT_DEFINITION_QUERY if not running_in_arc else None
    arcpy.management.MakeFeatureLayer(parks, "parks_lyr", definition_query)
    
    # Select and copy parks that intersect buffer
    arcpy.management.SelectLayerByLocation("parks_lyr", "INTERSECT", buffer_fc)
    arcpy.management.CopyFeatures("parks_lyr", out_fc)

    # Report results
    count = int(arcpy.management.GetCount(out_fc)[0])
    arcpy.AddMessage(f"Found {count} park polygons.")
    
    return out_fc


def get_parks_gdf_via_arc(point: Point, radius: float, parks: str = DEFAULT_WEB_LAYER, point_crs: str = "EPSG:4326") -> gpd.GeoDataFrame:
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
    fields = [f.name for f in arcpy.ListFields("parks_lyr") 
              if f.type not in ('Geometry', 'OID')]
    
    with arcpy.da.SearchCursor("parks_lyr", ["SHAPE@WKT"] + fields) as cursor:
        for row in cursor:
            geom = wkt.loads(row[0])
            geometries.append(geom)
            attributes.append(dict(zip(fields, row[1:])))
    
    # Create GeoDataFrame with the parks layer's CRS (not the input point's CRS)
    gdf = gpd.GeoDataFrame(attributes, geometry=geometries, crs=f"EPSG:{parks_epsg}")
    
    return gdf


def get_parks_gdf(point: Point, radius: float, parks: str = DEFAULT_WEB_LAYER, point_crs: str = "EPSG:3857") -> gpd.GeoDataFrame:
    """Get parks within a radius of a point.
    
    Args:
        point: Shapely Point geometry
        radius: Buffer radius in meters
        parks: URL to parks feature service
        point_crs: CRS of the input point (default: EPSG:3857)
    """
    
    # Convert point to WGS84
    point_wgs84 = gpd.GeoSeries([point], crs=point_crs).to_crs("EPSG:4326")[0]
    
    # Convert radius from meters to approximate degrees (~111km per degree)
    degree_radius = radius / 111000
    
    # Build bbox in degrees for server-side filter
    bbox = point_wgs84.buffer(degree_radius * 2).bounds
    
    # Server query (everything in WGS84)
    query_url = (
        f"{parks}/query?where={quote(DEFAULT_DEFINITION_QUERY)}"
        f"&geometry={bbox[0]},{bbox[1]},{bbox[2]},{bbox[3]}"
        f"&geometryType=esriGeometryEnvelope&inSR=4326"
        f"&spatialRel=esriSpatialRelIntersects"
        f"&outFields=*&f=geojson"
    )
    
    gdf = gpd.read_file(query_url)  # Already in WGS84
    
    if len(gdf) == 0:
        print("No parks found in query area")
        return gdf
    
    # Filter with degree-based buffer
    parks_in_radius = gdf[gdf.intersects(point_wgs84.buffer(degree_radius))]
    
    print(f"Found {len(parks_in_radius)} park polygons.")
    
    return parks_in_radius

if __name__ == "__main__":
    # Try to get ArcGIS tool parameters
    point_param = arcpy.GetParameter(0)
    
    # Check if point_param is a valid ArcGIS feature set (not None and has features)
    try:
        if point_param is not None:
            # Try to count features - if this works, we're in ArcGIS mode
            count = int(arcpy.management.GetCount(point_param)[0])
            if count > 0:
                # Running as ArcGIS tool with valid point
                point = point_param
                radius = arcpy.GetParameter(1)
                parks = arcpy.GetParameter(2)
                out_fc = arcpy.GetParameterAsText(3)
                arcpy.AddMessage(f"Running as ArcGIS tool with point: {point}, radius: {radius}, parks: {parks}, out_fc: {out_fc}")
            else:
                raise ValueError("No features in point parameter")
        else:
            raise ValueError("No point parameter")
    except Exception:
        # Running as standalone CLI
        point, radius, parks, out_fc = parse_cli_args()

    get_parks(point, radius, parks, out_fc)
