def interpolate_raster_nearest(dtm_path, dsm_path, lvm_path, mvm_path, hvm_path, bm_path, wm_path, bdm_path, 
                               igm_path, sm_path, tem_path, pnc_path, output_path):
    """
    Performs nearest-neighbor interpolation on a Digital Terrain Model (DTM) raster, filling in missing data points, and generates a classification map.

    **Steps:**

    1. **Loading and Reprojecting Rasters:**
       * Opens the DTM, DSM, and various classification masks (lvm, mvm, etc.).
       * Ensures all masks are reprojected to match the DTM's spatial reference and resolution.

    2. **Creating 'No-Data' Masks:**
       * Generates binary masks for each raster where 0 represents 'no-data' and a positive integer represents pixels of that specific classification.

    3. **Combining 'No-Data' Masks:**
       * Creates a composite 'no-data' mask. Pixels with a value of 0 in this mask are truly 'no-data' (not classified by any of the input masks).

    4. **DTM Interpolation:**
       * Interpolates the DTM using a nearest-neighbor algorithm, replacing 'no-data' values with the nearest valid neighboring pixel.
       * Converts the interpolated DTM to meters using a conversion factor.

    5. **Writing Interpolated DTM:**
        * Re-introduces 'no-data' (NaN) values in the interpolated DTM based on the combined 'no-data' mask.
        * Saves the final interpolated DTM to the `output_path`.

    6. **Generating Classification Map:**
       * **Identifies pixels as either 'unclassified' or 'ground' using a threshold comparison between DTM and DSM.**  
       * **Determines maximum vegetation or building class:** Calculates the maximum difference between DSM and various vegetation/building masks and assigns corresponding classification labels.
       * **Creates Classification Raster:** Constructs a multi-band raster representing the different 'no-data' classes, allowing for further analysis. 
       * **Saves Classification Raster:** Writes the classification map as a GeoTIFF file.

    7. **Cleanup:**
       * Closes raster datasets (if necessary for the format).
       * Deletes variables and calls the garbage collector for memory management.
    """

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