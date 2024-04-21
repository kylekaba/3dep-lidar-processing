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