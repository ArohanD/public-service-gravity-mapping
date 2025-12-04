import arcpy
from shapely import Polygon, wkt

def update_map(out_fc):
    """Update the map with the new demand feature class."""
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        m = aprx.activeMap
        m.addDataFromPath(out_fc)
        arcpy.AddMessage("Updated map with new demand feature class.")
    except Exception as e:
        arcpy.AddError(f"Error updating map: {e}")
        pass

def load_first_polygon_from_shapefile(shapefile: str) -> Polygon:
    """Load the first polygon from a shapefile into a polygon."""
    with arcpy.da.SearchCursor(shapefile, ["SHAPE@WKT"]) as cursor:
        row = next(cursor)
        polygon = wkt.loads(row[0])
    # verify that the polygon is valid using shapely
    if not polygon.is_valid:
        raise ValueError("Invalid polygon provided")
    return polygon