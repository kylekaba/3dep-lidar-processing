# lidar-utils: An aerial lidar point cloud data processing pipeline 
The programs within this directory were designed with the intent of processing aerial lidar data point clouds taken from the United States Geological Survey's 3D Elevation Program (3DEP). The code is designed to receive a path to an input directory that contains the .las or .laz point cloud files, a path to an output directory where the derivative products will be stored, the number of CPU cores that will be used, and keywords that can specify whether the entire directory or a subset of files will be processed.

## Install the Environment
First, install the necessary Python packages through the environment.yml file. To do this, download the repository and perform the following commands.

```
cd 3dep-lidar-processing
conda env create -f environment.yml
```
