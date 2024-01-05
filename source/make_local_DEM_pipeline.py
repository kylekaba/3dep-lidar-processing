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
    
    writer = {
        "type": "writers.gdal",
        "filename": output_path,
        "output_type": "idw" if dem_type == 'dsm' else "min",
        "resolution": dem_resolution
    }
    
    pipeline_stages = [reader] + filters + [writer]
    pipeline_dict = {"pipeline": pipeline_stages}
    
    return json.dumps(pipeline_dict)
