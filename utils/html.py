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
        .park { margin: 30px 0; padding: 20px; border: 1px solid #ddd; }
        .park h2 { margin-top: 0; }
        .park-content { display: flex; gap: 20px; }
        .park-map { flex: 1; height: 250px; }
        .park-stats { flex: 1; }
        table { width: 100%; border-collapse: collapse; }
        th, td { text-align: left; padding: 8px; border-bottom: 1px solid #ddd; }
        .negative { color: #c00; }
        .positive { color: #080; }
    </style>
</head>
<body>
    <h1>Park Demand Analysis Report</h1>
    <p>Top {{ parks|length }} parks most impacted by development (by weighted acres per 1000 change)</p>
    
    {% for park in parks %}
    <div class="park">
        <h2>#{{ park.rank }}: {{ park.name }}</h2>
        <div class="park-content">
            <div id="map-{{ park.rank }}" class="park-map"></div>
            <div class="park-stats">
                <table>
                    <tr>
                        <th>Metric</th>
                        <th>Current</th>
                        <th>Projected</th>
                        <th>Change</th>
                    </tr>
                    <tr>
                        <td>Population (weighted)</td>
                        <td>{{ "{:,.0f}".format(park.current_pop_weighted) }}</td>
                        <td>{{ "{:,.0f}".format(park.projected_pop_weighted) }}</td>
                        <td class="{% if park.pop_weighted_change > 0 %}negative{% else %}positive{% endif %}">
                            {{ "{:+,.0f}".format(park.pop_weighted_change) }}
                        </td>
                    </tr>
                    <tr>
                        <td>Acres per 1000 (weighted)</td>
                        <td>{{ "{:.2f}".format(park.current_acres_per_1000_weighted) }}</td>
                        <td>{{ "{:.2f}".format(park.projected_acres_per_1000_weighted) }}</td>
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
    """
    parks_gdf = parks_gdf.copy()
    parks_gdf["_change"] = parks_gdf["projected_acres_per_1000_weighted"].fillna(
        0
    ) - parks_gdf["current_acres_per_1000_weighted"].fillna(0)

    parks_with_change = parks_gdf[parks_gdf["_change"] != 0]
    top_parks = parks_with_change.nsmallest(top_n, "_change")

    park_reports = []
    for rank, (idx, row) in enumerate(top_parks.iterrows(), 1):
        park_reports.append(
            {
                "rank": rank,
                "name": _get_park_name(row),
                "current_pop_weighted": _safe_get(row, "current_pop_weighted"),
                "projected_pop_weighted": _safe_get(row, "projected_pop_weighted"),
                "pop_weighted_change": _safe_get(row, "pop_weighted_change"),
                "current_acres_per_1000_weighted": _safe_get(
                    row, "current_acres_per_1000_weighted"
                ),
                "projected_acres_per_1000_weighted": _safe_get(
                    row, "projected_acres_per_1000_weighted"
                ),
                "acres_weighted_change": _safe_get(
                    row, "acres_per_1000_weighted_change"
                ),
                "acres_weighted_pct_change": _safe_get(
                    row, "acres_per_1000_weighted_pct_change"
                ),
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
