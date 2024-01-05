import os
import json
import pdal
import rioxarray as rio
from pathlib import Path


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

        # Generate BVM (Building and Vegetation Mask)
        bvm_output_path = str(output_directory / f"{base_name}_bvm.tif")
        bvm_pipeline_json = make_local_DEM_pipeline(str(laz_file), 'bvm', 1.0, dsm_output_path)
        pdal.Pipeline(bvm_pipeline_json).execute()

        # Create CHM
        chm_output_path = str(output_directory / f"{base_name}_chm.tif")
        create_chm(dsm_output_path, dtm_output_path, chm_output_path)

        print(f"Processed DTM, DSM, BVM, and CHM for: {laz_file.name}")