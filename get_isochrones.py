import arcpy
import json
import urllib.request
import urllib.parse

from config import MAPBOX_TOKEN

def fetch_isochrone(lon, lat, token=MAPBOX_TOKEN) -> dict:
    """
    Fetch an isochrone from the Mapbox API for a given longitude and latitude.
    Returns a JSON object containing the isochrone polygons, or None on error.
    """
    base_url = "https://api.mapbox.com/isochrone/v1/mapbox/driving-traffic"
    coords = f"{lon},{lat}"

    params = {
        "contours_minutes": "10",
        "polygons": "true",
        "denoise": "1",
        "access_token": token
    }

    url = f"{base_url}/{coords}?{urllib.parse.urlencode(params)}"
    arcpy.AddMessage(f"Requesting isochrone: {url}")

    try:
        with urllib.request.urlopen(url) as response:
            text = response.read().decode("utf-8")
        parsed_response = json.loads(text)
        feature = parsed_response["features"][0]
        return feature
    except Exception as e:
        arcpy.AddError(f"Failed to fetch isochrone for ({lon}, {lat}): {e}")
        return None


# Script tool params
parks_layer = arcpy.GetParameter(0)        # Feature Layer of parks
token = arcpy.GetParameterAsText(1)        # Mapbox token
out_fc = arcpy.GetParameterAsText(2)       # Output feature class

# Defaults for standalone execution
if not parks_layer:
    parks_layer = r"..\Data\samples\OSM_NA_Leisure_GetParks_50.shp"

wgs84 = arcpy.SpatialReference(4326)
iso_geoms = []

def get_centroid_lonlat(geom, target_sr):
    """
    geom: Polygon geometry
    target_sr: SpatialReference to project to (WGS84)
    Returns (lon, lat)
    """
    # centroid is a Point (not a PointGeometry)
    centroid = geom.centroid    # Point

    # Wrap in PointGeometry so we can project
    point_geom = arcpy.PointGeometry(
        arcpy.Point(centroid.X, centroid.Y),
        geom.spatialReference
    )

    # Project the point geometry
    centroid_wgs = point_geom.projectAs(target_sr)  # PointGeometry

    # Get the underlying Point and return its X/Y
    pt = centroid_wgs.firstPoint    # Point
    return pt.X, pt.Y

if __name__ == "__main__":
    with arcpy.da.SearchCursor(parks_layer, ["SHAPE@"]) as cursor:
        for (park_geom,) in cursor:
            lon, lat = get_centroid_lonlat(park_geom, wgs84)
            arcpy.AddMessage(f"Park centroid: lon={lon}, lat={lat}")

            iso_data = fetch_isochrone(lon, lat)
            
            if iso_data is None:
                continue

            geojson_geom = iso_data["geometry"]
            g = arcpy.AsShape(geojson_geom)
            iso_geoms.append(g)

    if iso_geoms:
        arcpy.AddMessage(f"Created {len(iso_geoms)} isochrone polygons total.")
        
        if not out_fc:
            # Running standalone - save to scratch
            out_fc = r"C:\gispy\scratch\OSM_NA_Leisure_GetParks_50_getIsochrones.shp"
            arcpy.CopyFeatures_management(iso_geoms, out_fc)
            arcpy.AddMessage(f"Saved to: {out_fc}")
        else:
            # Running as a tool - save to specified output and add to map
            arcpy.CopyFeatures_management(iso_geoms, out_fc)
            aprx = arcpy.mp.ArcGISProject("CURRENT")
            m = aprx.activeMap
            m.addDataFromPath(out_fc)
    else:
        arcpy.AddWarning("No isochrones were created.")



