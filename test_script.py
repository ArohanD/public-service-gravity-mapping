import arcpy

feat = arcpy.GetParameter(0)

with arcpy.da.SearchCursor(feat, ["SHAPE@XY"]) as cursor:
    cursor = arcpy.da.SearchCursor(feat, ["SHAPE@XY"])
    row = next(cursor)
    x, y = row[0]

    print(x, y)
    arcpy.AddMessage(f"{x}, {y}")

# clear features so the layer is empty next time
arcpy.management.DeleteRows(feat)

# try to remove the layer from the active map
aprx = arcpy.mp.ArcGISProject("CURRENT")
m = aprx.activeMap
try:
    m.removeLayer(feat)
except:
    # If 'feat' is a pure Feature Set and not a layer, removeLayer will fail – that's OK.
    pass
