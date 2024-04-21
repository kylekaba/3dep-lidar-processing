def process_directory(directory_path, output_directory, num_cores, n_start=None, n_end=None):
    """
    Coordinates the processing of multiple LAZ files within a directory, designed for parallel execution.

    Args:
        directory_path (str): Path to the directory containing LAZ files.
        output_directory (str): Path where the results will be stored.
        num_cores (int): Number of CPU cores to utilize for parallel processing.
        n_start (int, optional): Index of the first LAZ file to process (for slicing).
        n_end (int, optional): Index of the last LAZ file to process (for slicing).
    """
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
        pool.starmap(process_single_file, zip(laz_files_to_process, [oput_dir] * len(laz_files_to_process)))
        pool.close()
        pool.join()
    
    gc.collect()
