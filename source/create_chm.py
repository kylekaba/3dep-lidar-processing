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