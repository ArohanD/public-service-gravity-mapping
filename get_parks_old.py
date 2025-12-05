# Input a point and a radius (in meters) and then output polygons
# representing all parks that touch the radius.

import arcpy
import sys
from get_point import get_point

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
        coords = sys.argv[1].split(",")
        x, y = float(coords[0]), float(coords[1])
    else:
        x, y = DEFAULT_COORDS

    radius = float(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_RADIUS
    parks = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_WEB_LAYER
    out_fc = sys.argv[4] if len(sys.argv) > 4 else DEFAULT_OUT_FC

    return (x, y), radius, parks, out_fc


def get_parks(point, radius, parks, out_fc):
    """Get parks within a radius of a point. Returns a feature class of parks that intersect the buffer."""

    running_in_arc = not isinstance(point, tuple)

    # Clean up existing layers
    if arcpy.Exists("parks_lyr"):
        arcpy.management.Delete("parks_lyr")
    if arcpy.Exists("in_memory\\click_buffer"):
        arcpy.management.Delete("in_memory\\click_buffer")

    # Extract coordinates from point (handle both feature and tuple)
    x, y = get_point(point) if running_in_arc else point

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


def get_parks_debug(point, radius, parks, out_fc):
    """Get parks within a radius of a point. Returns a feature class of parks that intersect the buffer."""

    running_in_arc = not isinstance(point, tuple)

    # DEBUG: Clean up existing layers
    if arcpy.Exists("parks_lyr"):
        arcpy.management.Delete("parks_lyr")
        arcpy.AddMessage("Deleted existing parks_lyr")
    if arcpy.Exists("in_memory\\click_buffer"):
        arcpy.management.Delete("in_memory\\click_buffer")

    # Extract coordinates from point (handle both feature and tuple)
    x, y = get_point(point) if running_in_arc else point
    arcpy.AddMessage(f"DEBUG: Point coordinates: {x}, {y}")

    # Get spatial reference from the INPUT point (not parks layer) to avoid CRS mismatch
    if running_in_arc:
        sr = arcpy.Describe(point).spatialReference
    else:
        sr = arcpy.Describe(parks).spatialReference
    arcpy.AddMessage(f"DEBUG: Point spatial reference: {sr.name}")

    # Create point geometry and buffer in the point's CRS
    pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)
    buffer_fc = "in_memory\\click_buffer"
    arcpy.analysis.Buffer(pt_geom, buffer_fc, f"{radius} Meters")
    arcpy.AddMessage(f"DEBUG: Buffer extent: {arcpy.Describe(buffer_fc).extent}")

    # Create parks layer with definition query if using URL
    definition_query = DEFAULT_DEFINITION_QUERY if not running_in_arc else None
    arcpy.AddMessage(f"DEBUG: Definition query: {definition_query}")
    arcpy.management.MakeFeatureLayer(parks, "parks_lyr", definition_query)
    arcpy.AddMessage(
        f"DEBUG: Parks in layer: {arcpy.management.GetCount('parks_lyr')[0]}"
    )

    # Select and copy parks that intersect buffer
    arcpy.management.SelectLayerByLocation("parks_lyr", "INTERSECT", buffer_fc)
    arcpy.AddMessage(
        f"DEBUG: Parks after selection: {arcpy.management.GetCount('parks_lyr')[0]}"
    )

    arcpy.management.CopyFeatures("parks_lyr", out_fc)

    # Report results
    count = int(arcpy.management.GetCount(out_fc)[0])
    arcpy.AddMessage(f"Found {count} park polygons.")

    return out_fc


# Main execution
point_param = arcpy.GetParameter(0)

if point_param is not None:
    # Running as ArcGIS tool
    point = point_param
    radius = arcpy.GetParameter(1)
    parks = arcpy.GetParameter(2)
    out_fc = arcpy.GetParameterAsText(3)
else:
    # Running as standalone CLI
    point, radius, parks, out_fc = parse_cli_args()

get_parks_debug(point, radius, parks, out_fc)
