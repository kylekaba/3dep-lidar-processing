def create_empty_mask(dtm_path, output_path):
    """
    Creates an empty raster mask with the same dimensions, spatial reference, and 
    metadata as an input DTM (Digital Terrain Model) raster, with all values set to zero. 

    Args:
        dtm_path (str): Path to the input DTM raster file.
        output_path (str): Path where the empty mask raster file will be saved.
    """
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