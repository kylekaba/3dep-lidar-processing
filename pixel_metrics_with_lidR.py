import os
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import StrVector
import multiprocessing
from pathlib import Path

def generate_pixel_metrics(laz_file, output_dir, pixel_size, elevation_units, r_lib_path):
    """
    Generates specified pixel metrics for a given .laz file using the lidR package in R and saves them as individual GeoTIFF files.

    Args:
        laz_file (str): Path to the input .laz file.
        output_dir (str): Directory where the output will be saved.
        pixel_size (int): Desired pixel size for metrics calculation.
        elevation_units (str): Units of the elevation data ('m', 'ft', 'km').
        r_lib_path (str): Path to the R library.
    """
    try:
        # Set the R library path
        robjects.r(f'.libPaths("{r_lib_path}")')

        # Import R packages
        base = importr('base')
        utils = importr('utils')
        lidr = importr('lidR')
        raster = importr('raster')
        sf = importr('sf')
        terra = importr('terra')
        e1071 = importr('e1071')  # for skewness and kurtosis

        # Convert Python strings to R strings
        laz_file_r = StrVector([laz_file])
        base_name = Path(laz_file).stem

        # Define the metrics and their corresponding filenames
        metrics = {
            "zmax": "max(Z)",
            "zmean": "mean(Z)",
            "zsd": "sd(Z)",
            "zskew": "skewness(Z)",
            "zkurtosis": "kurtosis(Z)",
            "zmin": "min(Z)",
            "zq10": "quantile(Z, probs = 0.10)",
            "zq25": "quantile(Z, probs = 0.25)",
            "zq50": "quantile(Z, probs = 0.50)",
            "zq75": "quantile(Z, probs = 0.75)",
            "zq90": "quantile(Z, probs = 0.90)",
            "zq95": "quantile(Z, probs = 0.95)"
        }

        # Define R script to process the .laz file and save each metric as a GeoTIFF
        r_script = f"""
        library(lidR)
        library(raster)
        library(sf)
        library(terra)
        library(e1071)  # for skewness and kurtosis

        las <- readLAS('{laz_file}')
        if (lidR::is.empty(las)) stop('Empty LAS file: {laz_file}')

        # Convert elevation to meters if necessary
        conversion_factor <- 1
        if ('{elevation_units}' == 'ft') {{
          conversion_factor <- 0.3048
        }} else if ('{elevation_units}' == 'km') {{
          conversion_factor <- 1000
        }}
        las@data$Z <- las@data$Z * conversion_factor

        # Normalize heights
        dtm <- rasterize_terrain(las, res={pixel_size}, algorithm=knnidw(k=6L, p=2, rmax={2*pixel_size}))
        las <- normalize_height(las, dtm)

        # Define custom function to calculate metrics
        custom_metrics <- function(x) {{
          zmax <- max(x)
          zmean <- mean(x)
          zsd <- sd(x)
          zskew <- skewness(x)
          zkurtosis <- kurtosis(x)
          zmin <- min(x)
          zq10 <- quantile(x, probs = 0.10)
          zq25 <- quantile(x, probs = 0.25)
          zq50 <- quantile(x, probs = 0.50)
          zq75 <- quantile(x, probs = 0.75)
          zq90 <- quantile(x, probs = 0.90)
          zq95 <- quantile(x, probs = 0.95)
          return(list(zmax=zmax, zmean=zmean, zsd=zsd, zskew=zskew, zkurtosis=zkurtosis, zmin=zmin, zq10=zq10, zq25=zq25, zq50=zq50, zq75=zq75, zq90=zq90, zq95=zq95))
        }}

        metrics <- pixel_metrics(las, ~custom_metrics(Z), res={pixel_size})

        # Get the bounding box coordinates directly
        xmin <- min(las@data$X)
        ymin <- min(las@data$Y)
        xmax <- max(las@data$X)
        ymax <- max(las@data$Y)

        # Get the CRS information from the LAS file
        crs_info <- st_crs(las)

        # Save each metric as an individual GeoTIFF file with correct CRS and extent
        for (metric in names(metrics)) {{
          raster_metrics <- raster(metrics[[metric]])
          crs(raster_metrics) <- crs_info$wkt
          extent(raster_metrics) <- extent(c(xmin, xmax, ymin, ymax))
          writeRaster(raster_metrics, filename=paste0('{output_dir}/{base_name}_', metric, '_{pixel_size}m.tif'), format='GTiff', overwrite=TRUE)
        }}
        """

        # Execute the R script
        robjects.r(r_script)
        print(f"Metrics GeoTIFF generated for {laz_file} at {pixel_size}m resolution")

    except Exception as e:
        print(f"Error processing {laz_file}: {str(e)}")

def process_laz_files(input_dir, output_dir, pixel_size, elevation_units, rs_lib_path, num_cores=2):
    """
    Process multiple .laz files in parallel to generate specified pixel metrics GeoTIFFs.

    Args:
        input_dir (str): Directory containing .laz files.
        output_dir (str): Directory where the output will be saved.
        pixel_size (int): Desired pixel size for metrics calculation.
        elevation_units (str): Units of the elevation data ('m', 'ft', 'km').
        r_lib_path (str): Path to the R library.
        num_cores (int): Number of CPU cores to use for parallel processing. Default is 2.
    """
    # Get the list of .laz files
    laz_files = sorted(Path(input_dir).rglob('*.laz'))

    # Prepare arguments for parallel processing
    args = [(str(file), output_dir, pixel_size, elevation_units, r_lib_path) for file in laz_files]

    # Use multiprocessing to process files in parallel
    with multiprocessing.Pool(num_cores) as pool:
        pool.starmap(generate_pixel_metrics, args)

if __name__ == "__main__":
    # Example usage
    input_directory = 'path/to/laz/files'
    output_directory = 'path/to/output/directory'
    pixel_size = 10  # Custom pixel size in m
    elevation_units = 'm'  # User-defined elevation units ('m', 'ft', 'km')
    r_lib_path = "r-library-path"  # Path to your R library

    # Use 2 CPU cores for processing
    process_laz_files(input_directory, output_directory, pixel_size, elevation_units, r_lib_path, num_cores=2)
