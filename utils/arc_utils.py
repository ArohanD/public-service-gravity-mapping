import arcpy
from shapely import Polygon, wkt
import os
import tempfile
import geopandas as gpd


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


def gdf_to_featureclass(gdf, output_path: str = None, layer_name: str = "gdf_layer"):
    """
    Convert a GeoDataFrame to a feature class.

    Args:
        gdf: GeoDataFrame to convert
        output_path: Optional path for the feature class.
                     If None, creates in-memory feature class.
        layer_name: Name for the layer (used for in-memory path)

    Returns:
        Path to the created feature class
    """

    if output_path is None:
        output_path = f"memory/{layer_name}"

    # Clean up if exists
    if arcpy.Exists(output_path):
        arcpy.management.Delete(output_path)

    # Write GeoDataFrame to temporary GeoJSON, then convert to feature class
    with tempfile.NamedTemporaryFile(suffix=".geojson", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Drop any geometry columns that aren't the main geometry
        # (Arc doesn't like multiple geometry columns or GeoSeries columns)
        cols_to_drop = []
        for col in gdf.columns:
            if col == "geometry":
                continue
            # Check if column is a GeoSeries
            if isinstance(gdf[col], gpd.GeoSeries):
                cols_to_drop.append(col)
            # Or check if it contains shapely geometries
            elif len(gdf) > 0 and hasattr(gdf[col].iloc[0], "geom_type"):
                cols_to_drop.append(col)

        if cols_to_drop:
            arcpy.AddMessage(f"Dropping geometry columns for export: {cols_to_drop}")

        gdf_export = gdf.drop(columns=cols_to_drop, errors="ignore")

        gdf_export.to_file(tmp_path, driver="GeoJSON")
        arcpy.AddMessage(f"Wrote temp GeoJSON to {tmp_path}")

        arcpy.conversion.JSONToFeatures(tmp_path, output_path)
        arcpy.AddMessage(f"Created feature class at {output_path}")

    except Exception as e:
        arcpy.AddError(f"Error converting GeoDataFrame: {e}")
        raise
    finally:
        # Clean up temp file
        if os.path.exists(tmp_path):
            os.remove(tmp_path)

    return output_path


def add_gdf_to_map(gdf, layer_name: str = "Analysis Results"):
    """
    Convert a GeoDataFrame to a feature class and add it to the current ArcGIS Pro map.

    Args:
        gdf: GeoDataFrame to add
        layer_name: Name for the layer in the map

    Returns:
        Path to the created feature class
    """
    arcpy.AddMessage(
        f"Converting GeoDataFrame with {len(gdf)} features to feature class..."
    )

    # Create in-memory feature class
    fc_path = gdf_to_featureclass(gdf, layer_name=layer_name.replace(" ", "_"))

    # Add to current map
    try:
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        m = aprx.activeMap
        if m is None:
            arcpy.AddWarning(
                "No active map found. Feature class created but not added to map."
            )
        else:
            m.addDataFromPath(fc_path)
            arcpy.AddMessage(f"Added '{layer_name}' to map: {m.name}")
    except Exception as e:
        arcpy.AddWarning(f"Could not add to map: {e}")
        arcpy.AddMessage(f"Feature class saved to: {fc_path}")

    return fc_path


def load_first_polygon_from_shapefile(shapefile: str) -> Polygon:
    """Load the first polygon from a shapefile into a polygon."""
    with arcpy.da.SearchCursor(shapefile, ["SHAPE@WKT"]) as cursor:
        row = next(cursor)
        polygon = wkt.loads(row[0])
    # verify that the polygon is valid using shapely
    if not polygon.is_valid:
        raise ValueError("Invalid polygon provided")
    return polygon
