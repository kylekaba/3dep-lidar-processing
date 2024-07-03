#!/bin/bash

# Directory containing the .laz files
input_dir="./input"
output_dir="./output"

# Create the output directory if it doesn't exist
mkdir -p $output_dir

# Loop through each .laz file and start a Docker container for each
for file in $input_dir/*.laz; do
  filename=$(basename "$file")
  export FILE_NAME=$filename
  # Run a Docker container for each .laz file
  docker-compose run --rm lidar_processing /input/$FILE_NAME /output 
done

# --jobs 10
