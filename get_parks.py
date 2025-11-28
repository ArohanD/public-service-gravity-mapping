# Input a point and a radius (in meters) and then output polygons
# representing all parks that touch the radius.

import arcpy
import sys
from get_point import get_point

arcpy.env.overwriteOutput = True

# Default parameters for CLI usage
DEFAULT_WEB_LAYER = "https://services6.arcgis.com/Do88DoK2xjTUCXd1/arcgis/rest/services/OSM_NA_Leisure/FeatureServer/0"
DEFAULT_DEFINITION_QUERY = "leisure = 'park'"
DEFAULT_RADIUS = 50.0
DEFAULT_COORDS = (-8784459.1687, 4294077.940800004)
DEFAULT_OUT_FC = "c:\\gispy\\scratch\\parks_in_radius.shp"

arcpy.AddMessage(f"sys.argv: {sys.argv}")

# Check if running as ArcGIS tool or standalone CLI
point_param = arcpy.GetParameter(0)

if point_param is not None:
    # Running as ArcGIS tool - use GetParameter
    point = point_param
    radius = arcpy.GetParameter(1)
    parks = arcpy.GetParameter(2)
    out_fc = arcpy.GetParameterAsText(3)
    arcpy.AddMessage("Running as ArcGIS tool")
else:
    # Running standalone CLI - parse command line: python get_parks.py [x,y] [radius] [parks_url] [out_fc]
    if len(sys.argv) > 1:
        coords = sys.argv[1].split(',')
        x = float(coords[0])
        y = float(coords[1])
    else:
        x, y = DEFAULT_COORDS
    
    radius = float(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_RADIUS
    parks = sys.argv[3] if len(sys.argv) > 3 else DEFAULT_WEB_LAYER
    out_fc = sys.argv[4] if len(sys.argv) > 4 else DEFAULT_OUT_FC
    
    point = (x, y)  # Store as tuple for standalone mode
    arcpy.AddMessage("Running standalone CLI")


def get_parks(point, radius):
    """
    Get parks within a radius of a point. Returns a feature class of parks that intersect the buffer.
    """
    
    # Get XY - handle both feature (ArcGIS tool) and tuple (standalone CLI)
    if isinstance(point, tuple):
        # Standalone CLI mode: point is already (x, y)
        x, y = point
    else:
        # ArcGIS tool mode: extract from feature
        x, y = get_point(point)
    
    # Make a point geometry in the same spatial ref as the parks layer
    sr = arcpy.Describe(parks).spatialReference
    pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)

    # Buffer the point by <radius> meters
    buffer_fc = "in_memory\\click_buffer"
    arcpy.analysis.Buffer(pt_geom, buffer_fc, f"{radius} Meters")

    # Select parks that intersect the buffer
    parks_lyr = "parks_lyr"
    
    # Apply definition query when using web layer URL (standalone CLI mode)
    if isinstance(parks, str) and parks.startswith("http"):
        arcpy.management.MakeFeatureLayer(parks, parks_lyr, DEFAULT_DEFINITION_QUERY)
    else:
        # ArcGIS tool mode - layer may already have definition query applied
        arcpy.management.MakeFeatureLayer(parks, parks_lyr)
    
    arcpy.management.SelectLayerByLocation(
        in_layer=parks_lyr,
        overlap_type="INTERSECT",
        select_features=buffer_fc
    )

    # Copy selected parks to the output feature class
    arcpy.management.CopyFeatures(parks_lyr, out_fc)

    # Report how many parks we found
    count = int(arcpy.management.GetCount(out_fc)[0])
    arcpy.AddMessage(f"Found {count} park polygons.")

    return out_fc


# Run when the script is executed as a tool
get_parks(point, radius)
