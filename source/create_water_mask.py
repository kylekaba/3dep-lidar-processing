import os
import json
import pdal
import rioxarray as rio
from pathlib import Path
import numpy as np

# Add the imports required for creating the mask
from osgeo import gdal, osr

# Function to create a rasterized water mask
def create_water_mask(laz_file_path, water_mask_path, raster_template_path):
    # Use PDAL to filter out water points from the point cloud
    pipeline_json = json.dumps({
        "pipeline": [
            {
                "type": "readers.las",
                "filename": laz_file_path
            },
            {
                "type": "filters.range",
                "limits": "Classification[9:9]"  # Assuming 9 is the classification code for water
            },
            {
                "type": "writers.gdal",
                "filename": water_mask_path,
                "gdaldriver": "GTiff",
                "output_type": "all",
                "resolution": 1.0
            }
        ]
    })
    
    pdal.Pipeline(pipeline_json).execute()
    
    # Ensure the mask aligns with the DTM/DSM resolution and extent
    with rio.open_rasterio(raster_template_path) as template_raster:
        with rio.open_rasterio(water_mask_path) as water_mask:
            water_mask_matched = water_mask.rio.reproject_match(template_raster)
            water_mask_matched.rio.to_raster(water_mask_path)
    
    # The mask now has the same spatial reference, resolution, and extent as the template raster