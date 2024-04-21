# lidar-utils: An aerial lidar point cloud data processing pipeline 
The programs within this directory were designed with the intent of processing aerial lidar data point clouds taken from the United States Geological Survey's 3D Elevation Program (3DEP). The code is designed to receive a path to an input directory that contains the .las or .laz point cloud files, a path to an output directory where the derivative products will be stored, the number of CPU cores that will be used, and keywords that can specify whether the entire directory or a subset of files will be processed.

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
    process_directory(data_directory, oput_dir,3,0,10)
```

**data_directory**: The directory where the .las or .laz aerial lidar point clouds are hosted.
**oput_dir**: The desired output directory where the four lidar derivative products (Digital Terrain Model (DTM), Digital Surface Model (DSM), Canopy Height Model (CHM), and Land Cover Classification Maps.
