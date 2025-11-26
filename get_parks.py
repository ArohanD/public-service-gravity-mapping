# Input a point and a radius (in meters) and then output polygons
# representing all parks that touch the radius.

import arcpy
from get_point import get_point

# --- Script tool parameters ---
point = arcpy.GetParameter(0)              # click point (Feature Set / layer)
radius = arcpy.GetParameter(1)             # radius in meters (Double)
parks = arcpy.GetParameter(2)              # parks polygon layer
out_fc = arcpy.GetParameterAsText(3)       # output feature class path (string)

arcpy.AddMessage(f"Radius: {radius} meters")


def get_parks(point, radius):
    # Get XY from the helper
    x, y = get_point(point)
    arcpy.AddMessage(f"Clicked point: {x}, {y}")

    # Make a point geometry in the same spatial ref as the parks layer
    sr = arcpy.Describe(parks).spatialReference
    pt_geom = arcpy.PointGeometry(arcpy.Point(x, y), sr)

    # Buffer the point by <radius> meters
    buffer_fc = "in_memory\\click_buffer"
    arcpy.analysis.Buffer(pt_geom, buffer_fc, f"{radius} Meters")

    # Select parks that intersect the buffer
    parks_lyr = "parks_lyr"
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
