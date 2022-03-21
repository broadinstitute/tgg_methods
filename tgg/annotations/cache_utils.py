import os
import pandas as pd
import re

CACHE_DIR = os.path.expanduser("~/.annotations")


def cache_data_table(get_table_func):
    """Decorator that caches the pandas DataFrame returned by the decorated function.
    It's intended for functions that take a relatively long time to retrieve some table over the network.
    Before calling the decorated function, the decorator checks whether result already exists in the
    cache dir (~/.annotations). If yes, it just reads the table from disk and returns it.
    If no, it calls the function and then saves the result table to ~/.annotations before returning it.
    """

    def wrapper(*args, **kwargs):
        # create cache dir
        if not os.path.isdir(CACHE_DIR):
            os.mkdir(CACHE_DIR)

        # check if cached file already exists
        function_name = get_table_func.__name__
        filename = re.sub("^get_", "", function_name) + ".tsv.gz"
        cache_file_path = os.path.join(CACHE_DIR, filename)

        if os.path.isfile(cache_file_path):
            return pd.read_table(cache_file_path)

        # call the underlying function
        df = get_table_func(*args, **kwargs)

        # save result to cache
        df.to_csv(cache_file_path, header=True, index=False, sep="\t")

        return df

    return wrapper
