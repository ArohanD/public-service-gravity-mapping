import arcpy
import json
import urllib.request
import urllib.parse

from config import MAPBOX_TOKEN


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
