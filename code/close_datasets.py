def close_datasets(*args):
    """
    Close any open rioxarray datasets to free up memory.
    """
    for arg in args:
        if isinstance(arg, rio.rioxarray.RasterArray):
            arg.close()