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

def make_local_DEM_pipeline(laz_file_path, dem_type='dtm', dem_resolution=1.0, output_path='output_dem.tif', filterNoise=True, elevation_units='m', epsg_code='3857'):
    """
    Constructs a PDAL pipeline configuration (JSON) for generating various types of Digital Elevation Models (DEMs) or masks from a LAZ file.

    Args:
        laz_file_path (str): Path to the input LAZ file.
        dem_type (str, optional):  The type of DEM or mask to generate. Options include: 'dtm' (bare earth), 'dsm' (all points), 'pnc' (processed but unclassified), 'lvm' (low vegetation), 'mvm' (medium vegetation), 'hvm' (high vegetation), 'bm' (building), 'wm' (water), 'bdm' (bridge deck), 'igm' (ignored ground), 'sm' (snow), 'tem' (temporal exclusion). Defaults to 'dtm'.
        dem_resolution (float, optional):  Desired resolution of the output DEM in meters. Defaults to 1.0.
        output_path (str, optional): Path where the generated DEM GeoTIFF will be saved. Defaults to 'output_dem.tif'.
        filterNoise (bool, optional): If True, applies filters to remove noise points at the beginning of the pipeline. Defaults to True.
        elevation_units (str, optional): Units of the elevation data ('m', 'ft', 'km'). Defaults to 'm'.
        epsg_code (str, optional): EPSG code for the desired output coordinate reference system. Defaults to '3857'.

    Returns:
        str:  A JSON string representing the configured PDAL pipeline.
    """
    
    reader = {
        "type": "readers.las",
        "filename": laz_file_path,
    }
    
    filters = []
    if elevation_units == 'ft':
        filters.append({"type": "filters.ferry", "dimensions": "Z=>Z_scaled"})
        filters.append({"type": "filters.assign", "value": "Z = Z_scaled * 0.3048"})
    elif elevation_units == 'km':
        filters.append({"type": "filters.ferry", "dimensions": "Z=>Z_scaled"})
        filters.append({"type": "filters.assign", "value": "Z = Z_scaled * 1000"})
    
    if filterNoise:
        filters.extend([{
            "type": "filters.range",
            "limits": "Classification![7:7]",
            },{"type": "filters.range",
            "limits": "Classification![18:18]"
        }])
    
    if dem_type == 'pnc':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[1:1]"
        })
    elif dem_type == 'dsm':
        pass
    elif dem_type == 'dtm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[2:2]"
        })
    elif dem_type == 'lvm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[3:3]"
        })
    elif dem_type == 'mvm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[4:4]"
        })
    elif dem_type == 'hvm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[5:5]"
        })
    elif dem_type == 'bm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[6:6]"
        })
    elif dem_type == 'wm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[9:9]"
        })
    elif dem_type == 'bdm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[17:17]"
        })
    elif dem_type == 'igm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[20:20]"
        })
    elif dem_type == 'sm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[21:21]"
        })
    elif dem_type == 'tem':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[22:22]"
        })

    reprojection_filter = {
        "type": "filters.reprojection",
        "out_srs": f"EPSG:{epsg_code}"
    }

    pipeline_stages = [reader] + filters + [reprojection_filter]

    writer = {
        "type": "writers.gdal",
        "filename": output_path,
        "output_type": "idw" if dem_type in ['dtm'] else "max",
        "radius": 2*float(dem_resolution),
        "data_type": "float32",
        "resolution": dem_resolution,
        "override_srs": f"EPSG:{epsg_code}"
    }
    
    pipeline_dict = {"pipeline": pipeline_stages + [writer]}
    
    return json.dumps(pipeline_dict)

def interpolate_raster_nearest(dtm_path, dsm_path, lvm_path, mvm_path, hvm_path, bm_path, wm_path, bdm_path, 
                               igm_path, sm_path, tem_path, pnc_path, output_path):
    """
    Performs nearest-neighbor interpolation on a Digital Terrain Model (DTM) raster, filling in missing data points, and generates a classification map.
    """
    dtm = rio.open_rasterio(dtm_path,masked=True)
    dsm = rio.open_rasterio(dsm_path, masked=True).rio.reproject_match(dtm)

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

    combined_no_data_mask = np.zeros_like(dtm_no_data_binary)

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
    
    dtm_interpolated = dtm.rio.interpolate_na(method="nearest")

    interpolated_dtm_data = dtm_interpolated.data[0]

    interpolated_dtm_data[np.where(combined_no_data_mask == 0)] = np.nan
    dtm_interpolated.data[0] = interpolated_dtm_data
    dtm_interpolated.rio.to_raster(output_path)

    combined_no_data_mask = np.zeros_like(dtm_no_data_binary)

    unclass_or_ground = (pnc_no_data_binary == 1) | (dtm_no_data_binary == 2)
    diff_pnc = np.abs(dsm.data[0] - pnc.data[0])
    diff_dtm = np.abs(dsm.data[0] - dtm.data[0])

    combined_no_data_mask[unclass_or_ground] = 1
    combined_no_data_mask[unclass_or_ground & (diff_dtm <= 0.5)] = 2

    choices = np.stack([lvm.data[0], mvm.data[0], hvm.data[0], bm.data[0]], axis=0)

    values = np.array([3, 4, 5, 6])

    max_indices = np.zeros_like(dsm.data[0], dtype=int)

    non_nan_mask = np.any(~np.isnan(choices), axis=0)

    max_indices[non_nan_mask] = np.nanargmax(choices[:, non_nan_mask], axis=0)

    valid_replacement_mask = non_nan_mask & (
        (lvm_no_data_binary != 0) | 
        (mvm_no_data_binary != 0) | 
        (hvm_no_data_binary != 0) | 
        (bm_no_data_binary != 0)
    )

    new_values = values[max_indices]

    combined_no_data_mask = np.where(
        valid_replacement_mask,
        new_values,
        combined_no_data_mask
    )

    no_data = np.zeros_like(dtm.values) + combined_no_data_mask

    no_data_mask = xr.DataArray(
        data=no_data,
        dims=dtm.dims,
        coords=dtm.coords,
        attrs=dtm.attrs
    )

    binary_masks = [pnc_no_data_binary, dtm_no_data_binary, lvm_no_data_binary, 
                    mvm_no_data_binary, hvm_no_data_binary, bm_no_data_binary]

    qualifying_masks = [combined_no_data_mask]
    for mask in binary_masks:
        if np.any(mask):
            qualifying_masks.append(mask.astype(np.uint8))

    stacked_masks = np.stack(qualifying_masks, axis=0) 

    mask_array = xr.DataArray(
        data=stacked_masks,
        dims=["band", "y", "x"],
        coords={
            "band": range(1, len(qualifying_masks) + 1), 
            "y": dtm.coords["y"],
            "x": dtm.coords["x"] 
        },
        attrs=dtm.attrs
    )
    mask_array.rio.write_crs(dtm.rio.crs,inplace=True)
    mask_array.rio.write_nodata(0,inplace=True)
    mask_array.rio.to_raster(dtm_path[0:-7] + "classification_maps.tif",dtype='uint8')

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

    del dtm, dsm, lvm, mvm, hvm, bm,wm, bdm, igm, sm, tem, pnc

    gc.collect()

def create_chm(dsm_path, interpolated_dtm_path, chm_path):
    dtm = rio.open_rasterio(interpolated_dtm_path, masked=True)
    dsm = rio.open_rasterio(dsm_path, masked=True)

    dsm.rio.to_raster(dsm_path)

    dsm.close()

    dsm = rio.open_rasterio(dsm_path, masked=True)

    dsm = dsm.rio.reproject_match(dtm)

    dsm = dsm.assign_coords({"x": dtm.x, "y": dtm.y})
    
    chm = dsm - dtm

    chm.rio.set_nodata(dtm.rio.nodata, inplace=True)
    chm.rio.to_raster(chm_path)

    dtm.close()
    dsm.close()
    chm.close()

    del dsm, dtm, chm

    gc.collect()

def create_empty_mask(dtm_path, output_path):
    dtm = rio.open_rasterio(dtm_path)

    empty_mask_data = np.zeros_like(dtm.values)

    empty_mask = xr.DataArray(
        data=empty_mask_data,
        dims=dtm.dims,
        coords=dtm.coords,
        attrs=dtm.attrs
    )

    empty_mask.rio.write_nodata(0, inplace=True)
    empty_mask.rio.to_raster(output_path)

def check_classification_points(laz_file_path, classification):
    pipeline_json = json.dumps({
    "pipeline": [
    {"type": "readers.las", "filename": laz_file_path},
    {"type": "filters.range", "limits": f"Classification[{classification}:{classification}]"}
]
    })

    pipeline = pdal.Pipeline(pipeline_json)
    pipeline.execute()

    points = pipeline.arrays[0]
    num_points = len(points)

    if num_points == 0:
        return False
    else:
        return True

def close_datasets(*args):
    for arg in args:
        if isinstance(arg, rio.rioxarray.RasterArray):
            arg.close()

def process_single_file(laz_file_path, output_directory, dem_resolution, elevation_units, epsg_code):
    try:
        base_name = Path(laz_file_path).stem

        dtm_output_path = os.path.join(output_directory, f"{base_name}_dtm.tif")
        dtm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'dtm', dem_resolution, dtm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
        pdal.Pipeline(dtm_pipeline_json).execute()

        dsm_output_path = os.path.join(output_directory, f"{base_name}_dsm.tif")
        dsm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'dsm', dem_resolution, dsm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
        pdal.Pipeline(dsm_pipeline_json).execute()

        pnc_output_path = output_directory + "/" + str(f"{base_name}_pnc.tif")
        check_pnc = check_classification_points(str(laz_file_path),1) 
        if check_pnc == False:
            create_empty_mask(dtm_output_path,pnc_output_path)
        elif check_pnc == True:
            pnc_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'pnc', dem_resolution, pnc_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(pnc_pipeline_json).execute()

        lvm_output_path = output_directory + "/" + str(f"{base_name}_lvm.tif")
        check_lv = check_classification_points(str(laz_file_path),3) 
        if check_lv == False:
            create_empty_mask(dtm_output_path,lvm_output_path)
        elif check_lv == True:
            lvm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'lvm', dem_resolution, lvm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(lvm_pipeline_json).execute()

        mvm_output_path = output_directory + "/" + str(f"{base_name}_mvm.tif")
        check_mv = check_classification_points(str(laz_file_path),4) 
        if check_mv == False:
            create_empty_mask(dtm_output_path,mvm_output_path)
        elif check_mv == True:
            mvm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'mvm', dem_resolution, mvm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(mvm_pipeline_json).execute()

        hvm_output_path = output_directory +"/" + str(f"{base_name}_hvm.tif")
        check_hv = check_classification_points(str(laz_file_path),5) 
        if check_hv == False:
            create_empty_mask(dtm_output_path,hvm_output_path)
        elif check_hv == True:
            hvm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'hvm', dem_resolution, hvm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(hvm_pipeline_json).execute()

        bm_output_path = output_directory + "/" + str(f"{base_name}_bm.tif")
        check_b = check_classification_points(str(laz_file_path),6) 
        if check_b == False:
            create_empty_mask(dtm_output_path, bm_output_path)
        elif check_b == True:
            bm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'bm', dem_resolution, bm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(bm_pipeline_json).execute()

        wm_output_path = output_directory + "/" + str(f"{base_name}_wm.tif")
        check_w = check_classification_points(str(laz_file_path), 9)
        if check_w == False:
            create_empty_mask(dtm_output_path, wm_output_path)
        elif check_w == True:
            wm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'wm', dem_resolution, wm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(wm_pipeline_json).execute()

        bdm_output_path = output_directory + "/" + str(f"{base_name}_bdm.tif")
        check_bd = check_classification_points(str(laz_file_path), 17)
        if check_bd == False:
            create_empty_mask(dtm_output_path, bdm_output_path)
        elif check_bd == True:
            bdm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'bdm', dem_resolution, bdm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(bdm_pipeline_json).execute()

        igm_output_path = output_directory + "/" + str(f"{base_name}_igm.tif")
        check_ig = check_classification_points(str(laz_file_path), 20)
        if check_ig == False:
            create_empty_mask(dtm_output_path, igm_output_path)
        elif check_ig == True:
            igm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'igm', dem_resolution, igm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(igm_pipeline_json).execute()

        sm_output_path = output_directory + "/" + str(f"{base_name}_sm.tif")
        check_s = check_classification_points(str(laz_file_path), 21)
        if check_s == False:
            create_empty_mask(dtm_output_path, sm_output_path)
        elif check_s == True:
            sm_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'sm', dem_resolution, sm_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(sm_pipeline_json).execute()

        tem_output_path = output_directory + "/" + str(f"{base_name}_tem.tif")
        check_te = check_classification_points(str(laz_file_path), 22)
        if check_te == False:
            create_empty_mask(dtm_output_path, tem_output_path)
        elif check_te == True:
            tem_pipeline_json = make_local_DEM_pipeline(str(laz_file_path), 'tem', dem_resolution, tem_output_path, filterNoise=True, elevation_units=elevation_units, epsg_code=epsg_code)
            pdal.Pipeline(tem_pipeline_json).execute()

        interpolated_dtm_output_path = output_directory + "/" + str(f"{base_name}_dtm_interpolated.tif")
        interpolate_raster_nearest(dtm_output_path, dsm_output_path, lvm_output_path, mvm_output_path, hvm_output_path, bm_output_path, wm_output_path,
                                   bdm_output_path, igm_output_path, sm_output_path, tem_output_path, pnc_output_path, interpolated_dtm_output_path)

        chm_output_path = output_directory + "/" + str(f"{base_name}_chm.tif")
        create_chm(dsm_output_path, interpolated_dtm_output_path, chm_output_path)

        outputs = [dtm_output_path, lvm_output_path, mvm_output_path, hvm_output_path, bm_output_path, wm_output_path, bdm_output_path, igm_output_path, sm_output_path, tem_output_path, pnc_output_path]
        for k in outputs:
            if os.path.exists(k):
                os.remove(k)
            else:
                print("The file does not exist")

        gc.collect()

        return str(laz_file_path), None
    except Exception as e:
        traceback.print_exc()
        return str(laz_file_path), e

def already_processed_check(laz_file, output_directory):
    base_name = Path(laz_file).stem
    chm_file = Path(output_directory) / f"{base_name}_chm.tif"
    return chm_file.exists()

def process_directory(directory_path, output_directory, num_cores, dem_resolution=1.0, elevation_units='m', epsg_code='3857', n_start=None, n_end=None):
    laz_files = sorted(Path(directory_path).rglob('*.laz'))

    if n_start is not None and n_end is not None:
        laz_files = laz_files[n_start:n_end] 

    laz_files_to_process = [f for f in laz_files if not already_processed_check(f, output_directory)]
 
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=num_cores, maxtasksperchild=1)
        pool.starmap(process_single_file, [(str(file), output_directory, dem_resolution, elevation_units, epsg_code) for file in laz_files_to_process])
        pool.close()
        pool.join()

    gc.collect()

# Example Usage: Using 3 CPUs, process the first 9 .laz files within the directory input at 1 meter resolution with EPSG Code 3857
data_directory = 'path/to/data'
output_directory = 'path/to/output/directory'
dem_resolution = 1  # User-defined DEM resolution
elevation_units = 'm'  # The units the elevation data are measured in (converted to meters) ('m', 'ft', 'km')
epsg_code = '3857'  # User-defined EPSG code

if __name__ == '__main__':
    process_directory(data_directory, output_directory, num_cores=3, dem_resolution=dem_resolution, elevation_units=elevation_units, epsg_code=epsg_code, n_start=0, n_end=10)
