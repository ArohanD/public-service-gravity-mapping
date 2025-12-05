"""HTML Report Generator for Park Demand Analysis."""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from jinja2 import Template
import geopandas as gpd
import io
import base64


REPORT_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Park Demand Analysis</title>
    <style>
        body { font-family: Arial, sans-serif; max-width: 900px; margin: 0 auto; padding: 20px; }
        h1 { border-bottom: 2px solid #333; padding-bottom: 10px; }
        .park { margin: 30px 0; padding: 20px; border: 1px solid #ddd; }
        .park h2 { margin-top: 0; }
        .park-content { display: flex; gap: 20px; }
        .park-map { flex: 1; }
        .park-map img { width: 100%; }
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
            <div class="park-map">
                {{ park.map_html | safe }}
            </div>
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
</body>
</html>
"""


def _create_park_map(park_row, parks_crs) -> str:
    """Create a base64-encoded PNG map image for a park."""
    park_geom = park_row["geometry"]
    isochrone_geom = park_row["isochrone_polygon"]
    
    park_gdf = gpd.GeoDataFrame(geometry=[park_geom], crs=parks_crs)
    isochrone_gdf = gpd.GeoDataFrame(geometry=[isochrone_geom], crs="EPSG:4326")
    park_gdf_plot = park_gdf.to_crs("EPSG:4326")
    
    fig, ax = plt.subplots(figsize=(4, 4))
    
    isochrone_gdf.plot(ax=ax, facecolor='lightblue', edgecolor='blue', alpha=0.3, linewidth=1)
    park_gdf_plot.plot(ax=ax, facecolor='green', edgecolor='darkgreen', alpha=0.6, linewidth=1)
    
    bounds = isochrone_gdf.total_bounds
    pad = (bounds[2] - bounds[0]) * 0.1
    ax.set_xlim(bounds[0] - pad, bounds[2] + pad)
    ax.set_ylim(bounds[1] - pad, bounds[3] - pad + (bounds[3] - bounds[1]) * 0.2)
    ax.set_axis_off()
    
    legend_elements = [
        Patch(facecolor='green', edgecolor='darkgreen', alpha=0.6, label='Park'),
        Patch(facecolor='lightblue', edgecolor='blue', alpha=0.3, label='Walk zone')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8)
    plt.tight_layout()
    
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
    buf.seek(0)
    img_base64 = base64.b64encode(buf.read()).decode('utf-8')
    plt.close(fig)
    
    return f'<img src="data:image/png;base64,{img_base64}" alt="Park map">'


def _get_park_name(park_row) -> str:
    """Extract park name from row, with fallback."""
    for field in ["name", "NAME", "park_name"]:
        if field in park_row.index and park_row[field]:
            return str(park_row[field])
    return f"Park #{park_row.name}"


def _safe_get(row, col, default=0):
    """Get value from row, returning default if None or missing."""
    val = row.get(col, default)
    return val if val is not None else default


def generate_park_report(parks_gdf: gpd.GeoDataFrame, output_path: str = "park_demand_report.html", top_n: int = 5) -> str:
    """
    Generate an HTML report for top N most impacted parks.
    
    Parks are sorted by weighted acres per 1000 change (most negative first).
    """
    # Calculate change for sorting
    parks_gdf = parks_gdf.copy()
    parks_gdf["_change"] = (
        parks_gdf["projected_acres_per_1000_weighted"].fillna(0) - 
        parks_gdf["current_acres_per_1000_weighted"].fillna(0)
    )
    
    # Get parks with actual change, sorted by most negative
    parks_with_change = parks_gdf[parks_gdf["_change"] != 0]
    top_parks = parks_with_change.nsmallest(top_n, "_change")
    
    # Build list of park dictionaries for template
    park_reports = []
    for rank, (idx, row) in enumerate(top_parks.iterrows(), 1):
        park_reports.append({
            "rank": rank,
            "name": _get_park_name(row),
            "current_pop_weighted": _safe_get(row, "current_pop_weighted"),
            "projected_pop_weighted": _safe_get(row, "projected_pop_weighted"),
            "pop_weighted_change": _safe_get(row, "pop_weighted_change"),
            "current_acres_per_1000_weighted": _safe_get(row, "current_acres_per_1000_weighted"),
            "projected_acres_per_1000_weighted": _safe_get(row, "projected_acres_per_1000_weighted"),
            "acres_weighted_change": _safe_get(row, "acres_per_1000_weighted_change"),
            "acres_weighted_pct_change": _safe_get(row, "acres_per_1000_weighted_pct_change"),
            "map_html": _create_park_map(row, parks_gdf.crs)
        })
    
    # Render and save
    html_content = Template(REPORT_TEMPLATE).render(parks=park_reports)
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html_content)
    
    print(f"Report saved to: {output_path}")
    return output_path
