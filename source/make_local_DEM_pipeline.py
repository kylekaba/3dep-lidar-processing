import os
import json
import pdal
import rioxarray as rio
from pathlib import Path

# Define the modified make_local_DEM_pipeline function
def make_local_DEM_pipeline(laz_file_path, dem_type='dtm', dem_resolution=1.0, output_path='output_dem.tif'):
    reader = {
        "type": "readers.las",
        "filename": laz_file_path
    }
    
    filters = []
    if dem_type == 'dtm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[2:2]"
        })
    # Bit of a tongue and cheek designation
    # The idea is to generate a building and vegetation model of just building and all vegetation points

    else if dem_type =='bvm' :
        filters.append({
            "type": "filters.range",
            "limits": "Classification[3:7]"
        })
    
    writer = {
        "type": "writers.gdal",
        "filename": output_path,
        "output_type": "idw" if dem_type == 'dsm' else "min",
        "resolution": dem_resolution
    }
    
    pipeline_stages = [reader] + filters + [writer]
    pipeline_dict = {"pipeline": pipeline_stages}
    
    return json.dumps(pipeline_dict)
