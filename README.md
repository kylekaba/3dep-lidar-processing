# lidar-utils: An aerial lidar point cloud data processing pipeline 
The programs within this directory were designed with the intent of processing aerial lidar data point clouds taken from the United States Geological Survey's 3D Elevation Program (3DEP). The code is designed to receive a path to an input directory that contains the .las or .laz point cloud files, a path to an output directory where the derivative products will be stored, the number of CPU cores that will be used, and keywords that can specify whether the entire directory or a subset of files will be processed. lidar-utils makes use of open-source software such as **PDAL** and **GDAL** as well as Python's in-built **multiprocessing** module to efficiently create rasterized elevation products.

## Install the Environment
First, install the necessary Python packages through the environment.yml file. To do this, download the repository and perform the following commands. The pipeline was developed with Python version 3.12.0.

```
cd 3dep-lidar-processing
conda env create -f environment.yml
```
## How to Run
As of this version, the most straight forward approach to apply this code is to replace the following variables within the **lidar_utils.py** script:

```
# Example Usage: Using 3 CPUs, process the first 9 .laz files within the directory input
data_directory = 'path-to-laz-files'
oput_dir = 'path-to-output-directory'
if __name__ == '__main__':
    process_directory(data_directory, oput_dir,num_cores,n_start,n_end)
```

*data_directory*: The directory where the .las or .laz aerial lidar point clouds are hosted.

*oput_dir*: The desired output directory where the four lidar derivative products (Digital Terrain Model (DTM), Digital Surface Model (DSM), Canopy Height Model (CHM), and Land Cover Classification Maps.

*num_cores*: The number of CPU processors that will be used. Each processor will independently process a .las or .laz file.

*n_start*: Starting index used to indicate where in the sorted list of .laz files the processing will begin. Default is **None**.

*n_end*: Ending index used to indicate where in the sorted list of .laz files the processing will end. Default is **None**.

Note that if *n_start* and *n_end* are both set to **None**, all of the .laz files within the directory will be processed.

## Description of Derivative Products

DTM: A digital terrain model that is 


## Additional Resources
Canopy Height Models, Digital Surface Models & Digital Elevation Models - Work With LiDAR Data in Python: https://www.earthdatascience.org/courses/use-data-open-source-python/data-stories/what-is-lidar-data/lidar-chm-dem-dsm/#:~:text=Digital%20Terrain%20Model%20(or%20DTM,of%20objects%20above%20the%20ground.

OpenTopography 3DEP Workflows: https://github.com/OpenTopography/OT_3DEP_Workflows

