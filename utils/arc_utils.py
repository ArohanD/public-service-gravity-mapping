# utils/arc_utils.py
#
# Author: Arohan Dutt
# Date: December 2025
#
# Purpose: Various ArcGIS utilities for the project.
#

import arcpy
import os
import tempfile
import geopandas as gpd


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
        # (Arc can get weird with multiple geometry columns)
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
