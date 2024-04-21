def already_processed_check(laz_file,output_directory):
    """
    Checks whether a .laz file has already been processed through the pipeline by checking for the presence of a CHM with the same extension.
    """
    # Extract the base name without the file extension
    base_name = Path(laz_file).stem
    # Check if the corresponding CHM File exists
    chm_file = Path(output_directory) / f"{base_name}_chm.tif"
    return chm_file.exists()
