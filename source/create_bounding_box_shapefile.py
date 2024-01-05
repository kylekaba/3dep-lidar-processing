import geopandas as gpd
from osgeo import gdal, osr
import os
from shapely.geometry import box

def create_bounding_box_shapefile(raster_path, shapefile_path):
    # Open the raster file
    ds = gdal.Open(raster_path)
    gt = ds.GetGeoTransform()

    # Coordinates of the bounding box
    x_min = gt[0]
    x_max = gt[0] + (gt[1] * ds.RasterXSize)
    y_min = gt[3] + (gt[5] * ds.RasterYSize)
    y_max = gt[3]

    # Create a polygon from the bounding box coordinates
    bounding_box_polygon = box(x_min, y_min, x_max, y_max)

    # Create a GeoDataFrame with the polygon
    gdf = gpd.GeoDataFrame([1], geometry=[bounding_box_polygon], columns=['dummy'])

    # Set the coordinate reference system (CRS) of the GeoDataFrame
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    gdf.crs = srs.ExportToProj4()

    # Write the GeoDataFrame to a new shapefile
    gdf.to_file(shapefile_path)

# Paths to input raster and output shapefile
#raster_path = '/Users/kkabasar/Desktop/usgs-3dep/OT_WebPortal_Products/rasters_NCAL/output_be.tif'  
#shapefile_path = '/Users/kkabasar/Desktop/usgs-3dep/data_workflow_test/shape-converts/howell-reservoir.shp'  

create_bounding_box_shapefile(raster_path, shapefile_path)