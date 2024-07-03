library(lidR)
library(raster)
library(sf)
library(terra)
library(e1071)  # for skewness and kurtosis

input_directory <- '/input'
output_directory <- '/output'
laz_file <- 'USGS_LPC_CA_SantaClaraCounty_2020_A20_10759750.laz'  # Replace with actual file name
elevation_units <- 'ft'  # Replace with actual units ('ft', 'km', or 'm')
pixel_size <- 100  # Replace with desired pixel size
base_name <- 'test'  # Replace with desired base name for output files

print(paste0(input_directory,'/', laz_file))
print(getwd())
las <- readLAS(paste0(input_directory,'/', laz_file))
if (lidR::is.empty(las)) stop(paste('Empty LAS file:', laz_file))
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
print("hi")
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
    writeRaster(raster_metrics, filename=paste0(output_directory,'/{base_name}_', metric, '_{pixel_size}m.tif'), format='GTiff', overwrite=TRUE)
}}