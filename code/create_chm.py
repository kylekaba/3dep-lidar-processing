def create_chm(dsm_path, interpolated_dtm_path, chm_path):
    """
    Generates a Canopy Height Model (CHM) raster by subtracting a Digital Terrain Model (DTM) from a Digital Surface Model (DSM).

    Steps:
        1. Loads the DSM and interpolated DTM rasters.
        2. Converts the DSM's units (if necessary) from US Survey Feet to meters.
        3. Reprojects the DSM to ensure alignment with the DTM in terms of resolution and spatial reference.
        4. Calculates the CHM by subtracting the DTM from the DSM (pixel-wise).
        5. Sets 'nodata' values in the CHM based on the DTM's 'nodata'.
        6. Saves the resulting CHM raster.
        7. Performs cleanup by closing datasets and releasing memory. 
    """
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