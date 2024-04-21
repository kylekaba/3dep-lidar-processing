def make_local_DEM_pipeline(laz_file_path, dem_type='dtm', dem_resolution=1.0, output_path='output_dem.tif', filterNoise=True):
    """
    Constructs a PDAL pipeline configuration (JSON) for generating various types of Digital Elevation Models (DEMs) or masks from a LAZ file.

    Args:
        laz_file_path (str): Path to the input LAZ file.
        dem_type (str, optional):  The type of DEM or mask to generate. Options include: 'dtm' (bare earth), 'dsm' (all points), 'pnc' (processed but unclassified), 'lvm' (low vegetation), 'mvm' (medium vegetation), 'hvm' (high vegetation), 'bm' (building), 'wm' (water), 'bdm' (bridge deck), 'igm' (ignored ground), 'sm' (snow), 'tem' (temporal exclusion). Defaults to 'dtm'.
        dem_resolution (float, optional):  Desired resolution of the output DEM in meters. Defaults to 1.0.
        output_path (str, optional): Path where the generated DEM GeoTIFF will be saved. Defaults to 'output_dem.tif'.
        filterNoise (bool, optional): If True, applies filters to remove noise points at the beginning of the pipeline. Defaults to True.

    Returns:
        str:  A JSON string representing the configured PDAL pipeline.
    """
    
    reader = {
        "type": "readers.las",
        "filename": laz_file_path,
    }
    
    # Remove points of noise at the beginning 

    filters = [{
            "type": "filters.range",
            "limits": "Classification![7:7]",
            },{"type": "filters.range",
            "limits": "Classification![18:18]"
        }]


    # Class 1 is Processed, but Not Classified ('pnc' stands for "processed, not classified")
    if dem_type == 'pnc':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[1:1]"
        })

        # Exclude High Noise Points from the DSM!
    elif dem_type == 'dsm':
        # No specific class filtering for the DSM since noise has already been removed
        pass

    # Class 2 is Bare Earth or Ground
    elif dem_type == 'dtm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[2:2]"
        })

    # Class 3 is Low Vegetation, typically below 1m ('lvm' stands for "low-vegetation mask")
    elif dem_type == 'lvm': 
        filters.append({
            "type": "filters.range",
            "limits": "Classification[3:3]"
        })

    # Class 4 is Medium Vegetation, typically bewteen 1m - 10m ('mvm' stands for "medium-vegetation mask")
    elif dem_type == 'mvm': 
        filters.append({
            "type": "filters.range",
            "limits": "Classification[4:4]"
        })

    # Class 5 is High Vegetation, typically >10 m ('hvm' stands for "high-vegetation mask")
    elif dem_type == 'hvm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[5:5]"
        })

    # Class 6 is Buildings ('bm' stands for "building mask")
    elif dem_type == 'bm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[6:6]"
        })    

        
    # Class 9 is Water ('wm' stands for "water mask")
    elif dem_type == 'wm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[9:9]"
        })    
    
    # Class 17 is Bridge Decks ('bdm' stands for "bridge deck mask")
    elif dem_type == 'bdm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[17:17]"
    
        })    


    # Class 20 is Ignored Ground (bare earth near breaklines) ('igm' stands for "ignored ground mask")
    elif dem_type == 'igm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[20:20]"
        })

    # Class 21 is Snow ('sm' stands for "snow mask")
    elif dem_type == 'sm':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[21:21]"
        })
    
    # Class 22 is Temporal Exclusion ('tem' stands for "temporal exclusion mask")
    elif dem_type == 'tem':
        filters.append({
            "type": "filters.range",
            "limits": "Classification[22:22]"
        })

    # Add a reprojection filter here to change the CRS to EPSG:3857
    # This should set the units to meters by default
    reprojection_filter = {
        "type": "filters.reprojection",
        "out_srs": "EPSG:3857"
    }

    pipeline_stages = [reader] + filters + [reprojection_filter]

    # Define the writer stage
    writer = {
        "type": "writers.gdal",

        "filename": output_path,
        # Bare Earth Should Use IDW, everything else + DSM should use "max"
         "output_type": "idw" if dem_type in ['dtm'] else "max",
        # Use IDW for DTM but Maximum point 2 raster method for everything else
        "radius": 2*float(dem_resolution),
        "data_type": "float32", # Ensure output is 32-bit floating point
        "resolution": dem_resolution,
        "override_srs": "EPSG:3857" # Ensure the writer uses EPSG 3857

    }
    
    pipeline_dict = {"pipeline": pipeline_stages + [writer]}
    
    return json.dumps(pipeline_dict)