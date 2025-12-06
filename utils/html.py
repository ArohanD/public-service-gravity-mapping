# utils/html.py
#
# Author: Arohan Dutt
# Date: December 2025
#
# Purpose: HTML report generator for the project. Using jinja for templating
#          because I saw it was in Arc Pro by default. Leaflet for mapping.
#

import json
from jinja2 import Template
import geopandas as gpd
from shapely.geometry import mapping

from utils.utils import log


# NRPA 2024 Agency Performance Review thresholds (acres per 1000 residents)
# Source: https://www.nrpa.org/siteassets/research/2024-agency-performance-review.pdf
# Page 11. Note that this should be adjusted for population size. THe below are
# accurate for Durham and Raleigh, but smaller jurisdictions like Chapel Hill
# should use different values from the table in the report. 
RATING_THRESHOLDS = {
    "excellent": 18.0,  # 75th percentile
    "good": 10.0,       # 50th percentile (median)
    "fair": 5.0,        # 25th percentile
    # Below 5.0 = poor
}

RATING_COLORS = {
    "excellent": "#2e7d32",  # Dark green
    "good": "#66bb6a",       # Light green
    "fair": "#ffa726",       # Orange
    "poor": "#ef5350",       # Red
}


def _get_rating(acres_per_1000):
    """Get rating label based on NRPA thresholds."""
    if acres_per_1000 is None:
        return "unknown"
    if acres_per_1000 >= RATING_THRESHOLDS["excellent"]:
        return "excellent"
    elif acres_per_1000 >= RATING_THRESHOLDS["good"]:
        return "good"
    elif acres_per_1000 >= RATING_THRESHOLDS["fair"]:
        return "fair"
    else:
        return "poor"


REPORT_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Park Demand Analysis</title>
    <!-- Leaflet CSS and JavaScript -->
    <link rel="stylesheet" href="https://unpkg.com/leaflet@1.9.4/dist/leaflet.css" />
    <script src="https://unpkg.com/leaflet@1.9.4/dist/leaflet.js"></script>
    <style>
        body { font-family: Arial, sans-serif; max-width: 900px; margin: 0 auto; padding: 20px; }
        h1 { border-bottom: 2px solid #333; padding-bottom: 10px; }
        .methodology { background: #f5f5f5; padding: 15px; margin-bottom: 20px; border-radius: 5px; font-size: 0.9em; }
        .methodology h3 { margin-top: 0; }
        .methodology table { font-size: 0.9em; margin: 10px 0; }
        .methodology td, .methodology th { padding: 4px 12px; }
        .park { margin: 30px 0; padding: 20px; border: 1px solid #ddd; }
        .park h2 { margin-top: 0; display: flex; align-items: center; justify-content: space-between; }
        .park-title { flex: 1; }
        .park-ratings { display: flex; align-items: center; gap: 6px; }
        .park-content { display: flex; gap: 20px; }
        .park-map { flex: 1; height: 250px; }
        .park-stats { flex: 1; }
        table { width: 100%; border-collapse: collapse; }
        th, td { text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }
        .negative { color: #c00; }
        .positive { color: #080; }
        .rating-badge { 
            display: inline-block; 
            padding: 1px 5px; 
            border-radius: 2px; 
            font-size: 0.65em; 
            font-weight: bold;
            color: white;
            text-transform: uppercase;
        }
        .rating-arrow { font-size: 0.9em; color: #666; }
        .rating-excellent { background-color: #2e7d32; }
        .rating-good { background-color: #66bb6a; }
        .rating-fair { background-color: #ffa726; }
        .rating-poor { background-color: #ef5350; }
        .rating-unknown { background-color: #9e9e9e; }
        .density-excellent { color: #2e7d32; font-weight: bold; }
        .density-good { color: #66bb6a; font-weight: bold; }
        .density-fair { color: #ffa726; font-weight: bold; }
        .density-poor { color: #ef5350; font-weight: bold; }
    </style>
</head>
<body>
    <h1>Park Demand Analysis Report</h1>
    
    <div class="methodology">
        <h3>Rating Methodology</h3>
        <p>Park density ratings are based on the <strong>NRPA 2024 Agency Performance Review</strong>, which reports 
        national benchmarks for acres of parkland per 1,000 residents:</p>
        <table>
            <tr>
                <th>Rating</th>
                <th>Acres/1000</th>
                <th>NRPA Percentile</th>
            </tr>
            <tr>
                <td><span class="rating-badge rating-excellent">Excellent</span></td>
                <td>≥ 18.0</td>
                <td>Above 75th percentile</td>
            </tr>
            <tr>
                <td><span class="rating-badge rating-good">Good</span></td>
                <td>10.0 – 18.0</td>
                <td>50th – 75th percentile</td>
            </tr>
            <tr>
                <td><span class="rating-badge rating-fair">Fair</span></td>
                <td>5.0 – 10.0</td>
                <td>25th – 50th percentile</td>
            </tr>
            <tr>
                <td><span class="rating-badge rating-poor">Poor</span></td>
                <td>< 5.0</td>
                <td>Below 25th percentile</td>
            </tr>
        </table>
        <p><em>Source: <a href="https://www.nrpa.org/publications-research/parkmetrics/">NRPA Park Metrics</a></em></p>
    </div>
    
    <p>Top {{ parks|length }} parks most impacted by development (by weighted acres per 1000 change)</p>
    
    {% for park in parks %}
    <div class="park">
        <h2>
            <span class="park-title">#{{ park.rank }}: {{ park.name }}</span>
            <span class="park-ratings">
                <span class="rating-badge rating-{{ park.current_rating }}">{{ park.current_rating }}</span>
                <span class="rating-arrow">→</span>
                <span class="rating-badge rating-{{ park.projected_rating }}">{{ park.projected_rating }}</span>
            </span>
        </h2>
        <div class="park-content">
            <div id="map-{{ park.rank }}" class="park-map"></div>
            <div class="park-stats">
                <p><strong>Park Area:</strong> {% if park.park_acres is not none %}{{ "{:.2f}".format(park.park_acres) }} acres{% else %}N/A{% endif %}</p>
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Current</th>
                        <th>Projected</th>
                        <th>Change</th>
                    </tr>
                    <tr>
                        <td>Utilization (weighted population)</td>
                        <td>{{ "{:,.0f}".format(park.current_pop_weighted) }}</td>
                        <td>{{ "{:,.0f}".format(park.projected_pop_weighted) }}</td>
                        <td class="{% if park.pop_weighted_change > 0 %}negative{% else %}positive{% endif %}">
                            {{ "{:+,.0f}".format(park.pop_weighted_change) }}
                        </td>
                    </tr>
                    <tr>
                        <td>Density (Acres per 1000)</td>
                        <td class="density-{{ park.current_rating }}">{{ "{:.2f}".format(park.current_acres_per_1000_weighted) }}</td>
                        <td class="density-{{ park.projected_rating }}">{{ "{:.2f}".format(park.projected_acres_per_1000_weighted) }}</td>
                        <td class="{% if park.acres_weighted_change < 0 %}negative{% else %}positive{% endif %}">
                            {{ "{:+.2f}".format(park.acres_weighted_change) }} ({{ "{:+.1f}".format(park.acres_weighted_pct_change) }}%)
                        </td>
                    </tr>
                </table>
            </div>
        </div>
    </div>
    {% endfor %}
    
    <script>
    {% for park in parks %}
    (function() {
        const map = L.map('map-{{ park.rank }}');
        L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
            attribution: '&copy; OpenStreetMap'
        }).addTo(map);
        
        const isochrone = L.geoJSON({{ park.isochrone_geojson | safe }}, {
            style: { color: 'blue', fillColor: 'lightblue', fillOpacity: 0.3, weight: 1 }
        }).addTo(map);
        
        L.geoJSON({{ park.park_geojson | safe }}, {
            style: { color: 'darkgreen', fillColor: 'green', fillOpacity: 0.6, weight: 1 }
        }).addTo(map);
        
        map.fitBounds(isochrone.getBounds().pad(0.1));
    })();
    {% endfor %}
    </script>
</body>
</html>
"""


def _geom_to_geojson(geom, source_crs, target_crs="EPSG:4326") -> str:
    """Convert a geometry to GeoJSON string for Leaflet mapping, reprojecting if needed."""
    gdf = gpd.GeoDataFrame(geometry=[geom], crs=source_crs)
    if str(source_crs) != target_crs:
        gdf = gdf.to_crs(target_crs)
    return json.dumps(mapping(gdf.geometry.iloc[0]))


def _get_park_name(park_row) -> str:
    """
    Extract park name from row, with fallback. I noticed that sometimes certain
    fields (like the name here) are not always standard in the data we got back
    from the OpenStreetMap API.
    """
    for field in ["name", "NAME", "park_name"]:
        if field in park_row.index and park_row[field]:
            return str(park_row[field])
    return f"Park #{park_row.name}"


def _safe_get(row, col, default=0):
    """Get value from row, returning default if None or missing."""
    val = row.get(col, default)
    return val if val is not None else default


def generate_park_report(
    parks_gdf: gpd.GeoDataFrame,
    output_path: str = "park_demand_report.html",
    top_n: int = 5,
) -> str:
    """
    Generate an HTML report for top N most impacted parks.

    Parks are sorted by weighted acres per 1000 change (most negative first).
    Ratings are based on NRPA 2024 Agency Performance Review thresholds.
    """
    parks_gdf = parks_gdf.copy()
    parks_gdf["_change"] = parks_gdf["projected_acres_per_1000_weighted"].fillna(
        0
    ) - parks_gdf["current_acres_per_1000_weighted"].fillna(0)

    parks_with_change = parks_gdf[parks_gdf["_change"] != 0]
    top_parks = parks_with_change.nsmallest(top_n, "_change")

    park_reports = []
    for rank, (idx, row) in enumerate(top_parks.iterrows(), 1):
        # Convert park area from m² to acres (1 acre = 4047 m²)
        park_area_m2 = row.get("current_park_area_m2")
        park_acres = park_area_m2 / 4047 if park_area_m2 is not None else None

        # Get density values for rating
        current_density = _safe_get(row, "current_acres_per_1000_weighted")
        projected_density = _safe_get(row, "projected_acres_per_1000_weighted")

        park_reports.append(
            {
                "rank": rank,
                "name": _get_park_name(row),
                "park_acres": park_acres,
                "current_pop_weighted": _safe_get(row, "current_pop_weighted"),
                "projected_pop_weighted": _safe_get(row, "projected_pop_weighted"),
                "pop_weighted_change": _safe_get(row, "pop_weighted_change"),
                "current_acres_per_1000_weighted": current_density,
                "projected_acres_per_1000_weighted": projected_density,
                "acres_weighted_change": _safe_get(
                    row, "acres_per_1000_weighted_change"
                ),
                "acres_weighted_pct_change": _safe_get(
                    row, "acres_per_1000_weighted_pct_change"
                ),
                "current_rating": _get_rating(current_density),
                "projected_rating": _get_rating(projected_density),
                "park_geojson": _geom_to_geojson(row["geometry"], parks_gdf.crs),
                "isochrone_geojson": _geom_to_geojson(
                    row["isochrone_polygon"], "EPSG:4326"
                ),
            }
        )

    html_content = Template(REPORT_TEMPLATE).render(parks=park_reports)

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    log(f"Report saved to: {output_path}")
    return output_path
