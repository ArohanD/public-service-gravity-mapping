import arcpy
import json
import urllib.request
import urllib.parse

def fetch_isochrone(lon, lat, token):
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
        return json.loads(text)
    except Exception as e:
        arcpy.AddError(f"Failed to fetch isochrone for ({lon}, {lat}): {e}")
        return None


# Script tool params
parks_layer = arcpy.GetParameter(0)        # Feature Layer of parks
token = arcpy.GetParameterAsText(1)        # Mapbox token
out_fc = arcpy.GetParameterAsText(2)       # Output feature class

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

with arcpy.da.SearchCursor(parks_layer, ["SHAPE@"]) as cursor:
    for (park_geom,) in cursor:
        lon, lat = get_centroid_lonlat(park_geom, wgs84)
        arcpy.AddMessage(f"Park centroid: lon={lon}, lat={lat}")

        iso_data = fetch_isochrone(lon, lat, token)
        
        if iso_data is None:
            continue

        for feat in iso_data.get("features", []):
            geojson_geom = feat["geometry"]
            g = arcpy.AsShape(geojson_geom)
            iso_geoms.append(g)

if iso_geoms:
    arcpy.CopyFeatures_management(iso_geoms, out_fc)
    arcpy.AddMessage(f"Created {len(iso_geoms)} isochrone polygons total.")

    aprx = arcpy.mp.ArcGISProject("CURRENT")
    m = aprx.activeMap
    m.addDataFromPath(out_fc)
else:
    arcpy.AddWarning("No isochrones were created.")



