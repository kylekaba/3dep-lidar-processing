import os
import gc
import json
import time 
import pdal
import rioxarray as rio
from scipy.spatial import cKDTree
from pathlib import Path
import numpy as np
import rioxarray as rio
from scipy.spatial import cKDTree
import xarray as xr
import traceback 
import multiprocessing
from multiprocessing import Pool, cpu_count

# Define the modified make_local_DEM_pipeline function
def make_local_DEM_pipeline(laz_file_path, dem_type='dtm', dem_resolution=1.0, output_path='output_dem.tif',filterNoise=True):
    reader = {
        "type": "readers.las",
        "filename": laz_file_path,
    }
    
    # Remove points of noise at the beginning 

    filters = [{
            "type": "filters.range",
            "limits": "Classification![7:7]",
            },{"type": "filters.range",
            "limits": "Classification![18:18]"
        }]


    # Class 1 is Processed, but Not Classified ('pnc' stands for "processed, not classified")
    if dem_type == 'pnc':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[1:1]"
        })

        # Exclude High Noise Points from the DSM!
    elif dem_type == 'dsm':
        # No specific class filtering for the DSM since noise has already been removed
        pass

    # Class 2 is Bare Earth or Ground
    elif dem_type == 'dtm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[2:2]"
        })

    # Class 3 is Low Vegetation, typically below 1m ('lvm' stands for "low-vegetation mask")
    elif dem_type == 'lvm': 
        filters.append({
            "type": "filters.range",
            "limits": "Classification[3:3]"
        })

    # Class 4 is Medium Vegetation, typically bewteen 1m - 10m ('mvm' stands for "medium-vegetation mask")
    elif dem_type == 'mvm': 
        filters.append({
            "type": "filters.range",
            "limits": "Classification[4:4]"
        })

    # Class 5 is High Vegetation, typically >10 m ('hvm' stands for "high-vegetation mask")
    elif dem_type == 'hvm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[5:5]"
        })

    # Class 6 is Buildings ('bm' stands for "building mask")
    elif dem_type == 'bm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[6:6]"
        })    

        
    # Class 9 is Water ('wm' stands for "water mask")
    elif dem_type == 'wm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[9:9]"
        })    
    
    # Class 17 is Bridge Decks ('bdm' stands for "bridge deck mask")
    elif dem_type == 'bdm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[17:17]"
    
        })    


    # Class 20 is Ignored Ground (bare earth near breaklines) ('igm' stands for "ignored ground mask")
    elif dem_type == 'igm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[20:20]"
        })

    # Class 21 is Snow ('sm' stands for "snow mask")
    elif dem_type == 'sm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[21:21]"
        })
    
    # Class 22 is Temporal Exclusion ('tem' stands for "temporal exclusion mask")
    elif dem_type == 'tem':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[22:22]"
        })

    # Add a reprojection filter here to change the CRS to EPSG:3857
    # This should set the units to meters by default
    reprojection_filter = {
        "type": "filters.reprojection",
        "out_srs": "EPSG:3857"
    }

    pipeline_stages = [reader] + filters + [reprojection_filter]

    # Define the writer stage
    writer = {
        "type": "writers.gdal",

        "filename": output_path,
        # Bare Earth Should Use IDW, everything else + DSM should use "max"
         "output_type": "idw" if dem_type in ['dtm'] else "max",
        # Use IDW for DTM but Maximum point 2 raster method for everything else
        "radius": 2*float(dem_resolution),
        "data_type": "float32", # Ensure output is 32-bit floating point
        "resolution": dem_resolution,
        "override_srs": "EPSG:3857" # Ensure the writer uses EPSG 3857

    }
    
    pipeline_dict = {"pipeline": pipeline_stages + [writer]}
    
    return json.dumps(pipeline_dict)


def interpolate_raster_nearest(dtm_path,dsm_path,lvm_path,mvm_path,hvm_path,bm_path,wm_path,bdm_path,
                               igm_path,sm_path,tem_path,pnc_path,output_path):

    
    # Conversion factor between US survey feet and meters
    conversion_factor = 1200/3937
    
    # Open the DTM raster
    dtm = rio.open_rasterio(dtm_path,masked=True)
    
    # Open the DSM raster
    dsm = rio.open_rasterio(dsm_path, masked=True).rio.reproject_match(dtm)

    # Reproject all masks to match DTM's extent and resolution
    lvm = rio.open_rasterio(lvm_path,masked=True).rio.reproject_match(dtm)
    mvm = rio.open_rasterio(mvm_path,masked=True).rio.reproject_match(dtm)
    hvm = rio.open_rasterio(hvm_path,masked=True).rio.reproject_match(dtm)
    bm  = rio.open_rasterio(bm_path,masked=True).rio.reproject_match(dtm)
    wm  = rio.open_rasterio(wm_path,masked=True).rio.reproject_match(dtm)
    bdm = rio.open_rasterio(bdm_path,masked=True).rio.reproject_match(dtm)
    igm = rio.open_rasterio(igm_path,masked=True).rio.reproject_match(dtm)
    sm = rio.open_rasterio(sm_path,masked=True).rio.reproject_match(dtm)
    tem = rio.open_rasterio(tem_path,masked=True).rio.reproject_match(dtm)
    pnc = rio.open_rasterio(pnc_path,masked=True).rio.reproject_match(dtm)

    # Create binary masks for 'no-data' values in each raster
    # For each mask, the positive integer represents pixels with that classification
    # Cells with a value of '0' are considered 'no data' (of that type)
    pnc_no_data_binary = (~np.isnan(pnc.data[0])).astype(int)*1 
    dtm_no_data_binary = (~np.isnan(dtm.data[0])).astype(int)*2 
    lvm_no_data_binary = (~np.isnan(lvm.data[0])).astype(int)*3
    mvm_no_data_binary = (~np.isnan(mvm.data[0])).astype(int)*4
    hvm_no_data_binary = (~np.isnan(hvm.data[0])).astype(int)*5
    bm_no_data_binary = (~np.isnan(bm.data[0])).astype(int)*6
    wm_no_data_binary = (~np.isnan(wm.data[0])).astype(int)*9
    bdm_no_data_binary = (~np.isnan(bdm.data[0])).astype(int)*17
    igm_no_data_binary = (~np.isnan(igm.data[0])).astype(int)*20
    sm_no_data_binary = (~np.isnan(sm.data[0])).astype(int)*21
    tem_no_data_binary = (~np.isnan(tem.data[0])).astype(int)*22

    # Combine all binary masks to identify "No-Data" pixels that are not classified as any of the above
    # Any pixel in this combined mask with a value of "0" is considered "No-Data"
    combined_no_data_mask = np.zeros_like(dtm_no_data_binary)

    # Change the values of the combined_no_data_mask to the corresponding class designation
    combined_no_data_mask[np.where(pnc_no_data_binary == 1)] = 1
    combined_no_data_mask[np.where(dtm_no_data_binary == 2)] = 2
    combined_no_data_mask[np.where(lvm_no_data_binary == 3)] = 3
    combined_no_data_mask[np.where(mvm_no_data_binary == 4)] = 4
    combined_no_data_mask[np.where(hvm_no_data_binary == 5)] = 5
    combined_no_data_mask[np.where(bm_no_data_binary == 6)] = 6
    combined_no_data_mask[np.where(wm_no_data_binary == 9)] = 9
    combined_no_data_mask[np.where(bdm_no_data_binary == 17)] = 17
    combined_no_data_mask[np.where(igm_no_data_binary == 20)] = 20
    combined_no_data_mask[np.where(sm_no_data_binary == 21)] = 21
    combined_no_data_mask[np.where(tem_no_data_binary == 22)] = 22
    
    # Interpolate the DTM using nearest neighbor, filling 'no-data' with nearest valid values
    dtm_interpolated = dtm.rio.interpolate_na(method="nearest")

    # Convert the interpolated DTM data to a numpy array 
    interpolated_dtm_data = dtm_interpolated.data[0]

    # Write out the interpolated raster
    # Reset pixels that have'no-data' in all masks back to NaN
    interpolated_dtm_data[np.where(combined_no_data_mask == 0)] = np.nan
    dtm_interpolated.data[0] = interpolated_dtm_data*conversion_factor # Convert to Meters
    dtm_interpolated.rio.to_raster(output_path)

    # Updated Logic: Construction of 'combined_no_data_mask'
    combined_no_data_mask = np.zeros_like(dtm_no_data_binary)

    # 1 & 2 : Unclassified & Ground 
    unclass_or_ground = (pnc_no_data_binary == 1) | (dtm_no_data_binary == 2)
    diff_pnc = np.abs(dsm.data[0] - pnc.data[0])
    diff_dtm = np.abs(dsm.data[0] - dtm.data[0])

    # Set a 50 cm threshold between the DSM and DTM to be considered "ground"
    combined_no_data_mask[unclass_or_ground & (diff_pnc < diff_dtm)] = 1
    combined_no_data_mask[unclass_or_ground & (diff_dtm <= 0.5)] = 2

    # Calculate absolute differences with DSM
    choices = np.stack([lvm.data[0], mvm.data[0], hvm.data[0], bm.data[0]], axis=0)

    # Define the corresponding class values
    values = np.array([3, 4, 5, 6])

    # Initialize an array to store indices of the maximum differences
    max_indices = np.zeros_like(dsm.data[0], dtype=int)

    # Mask to identify locations with at least one non-NaN value
    non_nan_mask = np.any(~np.isnan(choices), axis=0)

    # Calculate indices of the maximum non-NaN values only where it's possible
    max_indices[non_nan_mask] = np.nanargmax(choices[:, non_nan_mask], axis=0)

    # Mask for valid replacement: where there's at least one non-NaN and the mask is nonzero
    valid_replacement_mask = non_nan_mask & (
        (lvm_no_data_binary != 0) | 
        (mvm_no_data_binary != 0) | 
        (hvm_no_data_binary != 0) | 
        (bm_no_data_binary != 0)
    )

    # Determine new potential values based on the max diff indices
    new_values = values[max_indices]

    # Update combined_no_data_mask only where valid replacement mask is True
    combined_no_data_mask = np.where(
        valid_replacement_mask,
        new_values,
        combined_no_data_mask
    )

    # Create a new DataArray with the empty data and the same spatial coords as the DTM
    no_data = np.zeros_like(dtm.values) + combined_no_data_mask

    no_data_mask = xr.DataArray(
        data=no_data,
        dims=dtm.dims,
        coords=dtm.coords,
        attrs=dtm.attrs
    )

   # Binary masks for classes 1-6 
    binary_masks = [pnc_no_data_binary, dtm_no_data_binary, lvm_no_data_binary, 
                    mvm_no_data_binary, hvm_no_data_binary, bm_no_data_binary]

    # Filter out fully-zero masks
    qualifying_masks = [combined_no_data_mask]
    for mask in binary_masks:
        if np.any(mask):  # Check if any non-zero elements exist
            qualifying_masks.append(mask.astype(np.uint8))


    # Stack for multi-band output
    stacked_masks = np.stack(qualifying_masks, axis=0) 

    mask_array = xr.DataArray(
        data=stacked_masks,
        dims=["band", "y", "x"],  # Assuming standard y, x raster dimensions
        coords={
            "band": range(1, len(qualifying_masks) + 1), 
            "y": dtm.coords["y"],
            "x": dtm.coords["x"] 
        },
        attrs=dtm.attrs  # Carry over relevant attributes
)
    mask_array.rio.write_crs(dtm.rio.crs,inplace=True)
    mask_array.rio.write_nodata(0,inplace=True)
    mask_array.rio.to_raster(dtm_path[0:-7] + "classification_maps.tif",dtype='uint8')

    # Delete all the rasters
    # Close the datasets if they are in a format that requires it
    dsm.close()
    lvm.close()
    mvm.close()
    hvm.close()
    bm.close()
    wm.close()
    bdm.close()
    igm.close()
    sm.close()
    tem.close()
    pnc.close()

    # Clear the variables from memory
    del dtm, dsm, lvm, mvm, hvm, bm,wm, bdm, igm, sm, tem, pnc

    # Call garbage collector to free up memory
    gc.collect()

def create_chm(dsm_path, interpolated_dtm_path, chm_path):
    # Read DSM and DTM rasters
    dtm = rio.open_rasterio(interpolated_dtm_path, masked=True)
    dsm = rio.open_rasterio(dsm_path, masked=True)

    # Convert the DSM from units of US Survey Feet into Meters
    dsm.values *= 1200/3937

    dsm.rio.to_raster(dsm_path)

    # Close the DSM and Re-Open It
    dsm.close()

    dsm = rio.open_rasterio(dsm_path, masked=True)

    # Ensure that the rasters are on the same grid 
    dsm = dsm.rio.reproject_match(dtm)

    # Make sure that the extent of the DSM and DTM have the same coordinates
    dsm = dsm.assign_coords({"x": dtm.x, "y": dtm.y})
    
    # Perform raster subtraction to get CHM
    chm = dsm - dtm

    # Set the nodata value and write out the CHM raster
    chm.rio.set_nodata(dtm.rio.nodata, inplace=True)
    chm.rio.to_raster(chm_path)

    # Close the datasets if they are in a format that requires it
    dtm.close()
    dsm.close()
    chm.close()

    # Clear the variables from memory
    del dsm, dtm, chm

    # Call garbage collector to free up memory
    gc.collect()

def create_empty_mask(dtm_path, output_path):
    # Open the DTM raster
    dtm = rio.open_rasterio(dtm_path)


    # Create an empty mask with zeros, matching the DTM's shape and data model
    empty_mask_data = np.zeros_like(dtm.values)

    # Create a new DataArray with the empty data and the same spatial coords as the DTM
    empty_mask = xr.DataArray(
        data=empty_mask_data,
        dims=dtm.dims,
        coords=dtm.coords,
        attrs=dtm.attrs
    )

    # Set the nodata value to 0
    empty_mask.rio.write_nodata(0, inplace=True)
    # Write out the empty mask to a raster
    empty_mask.rio.to_raster(output_path)

# Function that checks whether a given classification has any points or not.
def check_classification_points(laz_file_path, classification):
    pipeline_json = json.dumps({
    "pipeline": [
    {"type": "readers.las", "filename": laz_file_path},
    {"type": "filters.range", "limits": f"Classification[{classification}:{classification}]"}
]
    })

    pipeline = pdal.Pipeline(pipeline_json)
    pipeline.execute()

    # Get the number of points after filtering
    points = pipeline.arrays[0]
    num_points = len(points)

    if num_points == 0:
        return False
    else:
        return True

def close_datasets(*args):
    """
    Close any open rioxarray datasets to free up memory.
    """
    for arg in args:
        if isinstance(arg, rio.rioxarray.RasterArray):
            arg.close()

def process_single_file(laz_file_path, output_directory):
  """
  Processes a single LiDAR LAZ file by generating DTM, DSM, CHM, and land classification map.

  Args:
      laz_file_path: Path to the input LAZ file.
      output_dir: Path to the output directory where results will be saved.
  """

  try:
    base_name = Path(laz_file_path).stem

    # Generate DTM (Digital Terrain Model)
    dtm_output_path = os.path.join(output_directory, f"{base_name}_dtm.tif")
    dtm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'dtm', 1.0, dtm_output_path, filterNoise=True)
    pdal.Pipeline(dtm_pipeline_json).execute()

    # Generate DSM (Digital Surface Model)
    dsm_output_path = os.path.join(output_directory, f"{base_name}_dsm.tif")
    dsm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'dsm', 1.0, dsm_output_path, filterNoise=True)
    pdal.Pipeline(dsm_pipeline_json).execute()

    # Generate PNC (Processed, but not classified Mask)
    pnc_output_path = output_directory + "/" + str(f"{base_name}_pnc.tif")
    #Check whether there are any points from Class 1 in the Point Cloud
    check_pnc = check_classification_points(str(laz_file_path),1) 
    if check_pnc == False:
        create_empty_mask(dtm_output_path,pnc_output_path)
    elif check_pnc == True:
        pnc_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'pnc', 1.0, pnc_output_path,filterNoise=True)
        pdal.Pipeline(pnc_pipeline_json).execute()

    # Generate LVM (Low Vegetation Mask)
    lvm_output_path = output_directory + "/" + str(f"{base_name}_lvm.tif")
    #Check whether there are any points from Class 3 in the Point Cloud
    check_lv = check_classification_points(str(laz_file_path),3) 
    if check_lv == False:
        create_empty_mask(dtm_output_path,lvm_output_path)
    elif check_lv == True:
        lvm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'lvm', 1.0, lvm_output_path,filterNoise=True)
        pdal.Pipeline(lvm_pipeline_json).execute()

    # Generate MVM (Medium Vegetation Mask)
    mvm_output_path = output_directory + "/" + str(f"{base_name}_mvm.tif")

    # Check whether there are any points from Class 4 in the Point Cloud
    check_mv = check_classification_points(str(laz_file_path),4) 
    if check_mv == False:
        create_empty_mask(dtm_output_path,mvm_output_path)
    elif check_mv == True:
        mvm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'mvm', 1.0, mvm_output_path,filterNoise=True)
        pdal.Pipeline(mvm_pipeline_json).execute()

    # Generate HVM (High Vegetation Mask)
    hvm_output_path = output_directory +"/" + str(f"{base_name}_hvm.tif")

    #Check whether there are any points from Class 5 in the Point Cloud
    check_hv = check_classification_points(str(laz_file_path),5) 
    if check_hv == False:
        create_empty_mask(dtm_output_path,hvm_output_path)
    elif check_hv == True:
        hvm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'hvm', 1.0, hvm_output_path,filterNoise=True)
        pdal.Pipeline(hvm_pipeline_json).execute()

    # Generate BM (Building Mask)
    bm_output_path = output_directory + "/" + str(f"{base_name}_bm.tif")

    #Check whether there are any points from Class 6 in the Point Cloud
    check_b = check_classification_points(str(laz_file_path),6) 
    if check_b == False:
        create_empty_mask(dtm_output_path,bm_output_path)
    elif check_b == True:
        bm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'bm', 1.0, bm_output_path,filterNoise=True)
        pdal.Pipeline(bm_pipeline_json).execute()

    # Generate WM (Water Mask)
    wm_output_path = output_directory + "/" + str(f"{base_name}_wm.tif")

    # Check whether there are any points from Class 9 in the Point Cloud
    check_w = check_classification_points(str(laz_file_path),9) 
    if check_w == False:
        create_empty_mask(dtm_output_path,wm_output_path)
    elif check_w == True:
        wm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'wm', 1.0, wm_output_path,filterNoise=True)
        pdal.Pipeline(wm_pipeline_json).execute()


    # Generate BDM (Bridge Decks Mask)
    bdm_output_path = output_directory + "/" + str(f"{base_name}_bdm.tif")

    # Check whether there are any points from Class 17 in the Point Cloud
    check_bd = check_classification_points(str(laz_file_path),17) 
    if check_bd == False:
        create_empty_mask(dtm_output_path,bdm_output_path)
    elif check_bd == True:
        bdm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'bdm', 1.0, bdm_output_path,filterNoise=True)
        pdal.Pipeline(bdm_pipeline_json).execute()

    # Generate IGM (Ignored Ground Mask)
    igm_output_path = output_directory + "/" + str(f"{base_name}_igm.tif")

    # Check whether there are any points from Class 20 in the Point Cloud
    check_ig = check_classification_points(str(laz_file_path),20) 
    if check_ig == False:
        create_empty_mask(dtm_output_path,igm_output_path)
    elif check_ig == True:
        igm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'igm', 1.0, igm_output_path,filterNoise=True)
        pdal.Pipeline(igm_pipeline_json).execute()

    # Generate SM (Snow Mask)
    sm_output_path = output_directory + "/" + str(f"{base_name}_sm.tif")

    # Check whether there are any points from Class 21 in the Point Cloud
    check_s = check_classification_points(str(laz_file_path),21) 
    if check_s == False:
        create_empty_mask(dtm_output_path,sm_output_path)
    elif check_s == True:
        sm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'sm', 1.0, sm_output_path,filterNoise=True)
        pdal.Pipeline(sm_pipeline_json).execute()

    # Generate TEM (Temporal Exclusion Mask)
    tem_output_path = output_directory + "/" + str(f"{base_name}_tem.tif")

    # Check whether there are any points from Class 22 in the Point Cloud
    check_te = check_classification_points(str(laz_file_path),22) 
    if check_te == False:
        create_empty_mask(dtm_output_path,tem_output_path)
    elif check_te == True:
        tem_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'tem', 1.0, tem_output_path,filterNoise=True)
        pdal.Pipeline(tem_pipeline_json).execute()

    # Interpolated DTM
    interpolated_dtm_output_path = output_directory + "/" + str(f"{base_name}_dtm_interpolated.tif")
    interpolate_raster_nearest(dtm_output_path,dsm_output_path,lvm_output_path,mvm_output_path,hvm_output_path,bm_output_path,wm_output_path,
                            bdm_output_path,igm_output_path,sm_output_path,tem_output_path,pnc_output_path,interpolated_dtm_output_path)

    # Create CHM using the interpolated DTM
    chm_output_path = output_directory + "/" + str(f"{base_name}_chm.tif")
    create_chm(dsm_output_path, interpolated_dtm_output_path, chm_output_path)

    # Delete every file except the .dsm, .chm, combined_mask, and interpolated dtm
    outputs = [dtm_output_path,lvm_output_path,mvm_output_path,hvm_output_path,bm_output_path,wm_output_path,bdm_output_path,igm_output_path,sm_output_path,tem_output_path,pnc_output_path]
    for k in outputs:
        if os.path.exists(k):
            os.remove(k)
        else:
            print("The file does not exist")

    # Call garbage collector
    gc.collect()

    return str(laz_file_path), None  # Successfully processed file
  except Exception as e:
      traceback.print_exc()
      return str(laz_file_path), e  # Error encountered

def already_processed_check(laz_file,output_directory):
    # Extract the base name without the file extension
    base_name = Path(laz_file).stem
    # Check if the corresponding CHM File exists
    chm_file = Path(output_directory) / f"{base_name}_chm.tif"
    return chm_file.exists()


def process_directory(directory_path,output_directory,num_cores,n_start=None,n_end=None):

    # Sort the laz files
    laz_files = sorted(Path(directory_path).rglob('*.laz'))

    # Apply slicing based on n_start and n_end
    if n_start is not None and n_end is not None:
        laz_files = laz_files[n_start:n_end] 
    
    # Filter out already processed files
    laz_files_to_process = [f for f in laz_files if not already_processed_check(f,output_directory)]
    
    # Use multiprocessing to map the pipeline to different files 
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=num_cores,maxtasksperchild=1)
        pool.starmap(process_single_file, zip(laz_files_to_process, [output_directory] * len(laz_files_to_process)))
        pool.close()
        pool.join()
    
    gc.collect()
