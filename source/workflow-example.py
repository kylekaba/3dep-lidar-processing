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
    
    writer = {
        "type": "writers.gdal",
        "filename": output_path,
        "output_type": "idw" if dem_type == 'dsm' else "min",
        "resolution": dem_resolution
    }
    
    pipeline_stages = [reader] + filters + [writer]
    pipeline_dict = {"pipeline": pipeline_stages}
    
    return json.dumps(pipeline_dict)

# Function to create a CHM from DSM and DTM
def create_chm(dsm_path, dtm_path, chm_path):
# Ensure that the dsm and dtm variables are not in memory
    try: 
        dsm.close()
        dtm.close()
        chm.close
        del chm
        del dsm
        del dtm
    except:
        pass

    # Read DSM and DTM rasters
    dsm = rio.open_rasterio(dsm_path, masked=True)
    dtm = rio.open_rasterio(dtm_path, masked=True)

	# Ensure that the rasters are the same shape
    if dsm.shape == dtm.shape:
	    print('dtm and dsm have the same shape. proceed to the next step.')
    elif dsm.shape > dtm.shape:
	    small = dtm_path
	    large = dsm_path
	    dsm = dsm.rio.reproject_match(dtm)
	    print("the smaller raster is {}. It has been used to clip the larger raster {}. proceed to the next step.".format(small,large))
    else:
	    small = dsm_path
	    large = dtm_path
	    dtm = dtm.rio.reproject_match(dsm)
	    print("the smaller raster is {}. It has been used to clip the larger raster {}. proceed to the next step.".format(small,large))

	# Make sure that the extent of the rasters have the same coordinates
    dsm = dsm.assign_coords({
    "x": dtm.x,
    "y": dtm.y,})
	
    # Perform raster subtraction to get CHM
    chm = dsm - dtm
    chm.compute()
    chm.rio.set_nodata(dtm.rio.nodata, inplace=True)
    chm.rio.to_raster(chm_path)

# Loop over all .laz files in the directory and create DSMs, DTMs, and CHMs
def process_directory(directory_path):
    for laz_file in Path(directory_path).rglob('*.laz'):
        base_name = laz_file.stem
        output_directory = laz_file.parent

        # Generate DTM
        dtm_output_path = str(output_directory / f"{base_name}_dtm.tif")
        dtm_pipeline_json = make_local_DEM_pipeline(str(laz_file), 'dtm', 1.0, dtm_output_path)
        pdal.Pipeline(dtm_pipeline_json).execute()

        # Generate DSM
        dsm_output_path = str(output_directory / f"{base_name}_dsm.tif")
        dsm_pipeline_json = make_local_DEM_pipeline(str(laz_file), 'dsm', 1.0, dsm_output_path)
        pdal.Pipeline(dsm_pipeline_json).execute()

        # Create CHM
        chm_output_path = str(output_directory / f"{base_name}_chm.tif")
        create_chm(dsm_output_path, dtm_output_path, chm_output_path)

        print(f"Processed DTM, DSM, and CHM for: {laz_file.name}")