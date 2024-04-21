# Function that checks whether a given classification has any points or not.
def check_classification_points(laz_file_path, classification):
    pipeline_json = json.dumps({
    "pipeline": [
    {"type": "readers.las", "filename": laz_file_path},
    {"type": "filters.range", "limits": f"Classification[{classification}:{classification}]"}
]
    })

    pipeline = pdal.Pipeline(pipeline_json)
    pipeline.execute()

    # Get the number of points after filtering
    points = pipeline.arrays[0]
    num_points = len(points)

    if num_points == 0:
        return False
    else:
        return True
