"""
Helper functions for all input/output related things.
"""

import os
import tarfile
import pickle
import json
import gzip
import shutil
import requests
from tqdm import tqdm

def save(obj, path):
    """ Saves an object to either pickle, json, or json.gz (determined by the extension in the file name).

    Parameters
    ----------
    obj:
        The object to be saved.
    path:
        The path to save the object.
    """
    if path.endswith('.json.gz'):
        with gzip.open(path, 'w') as file:
            file.write(json.dumps(obj).encode('utf-8'))
    elif path.endswith('.json'):
        with open(path,'w') as file:
            json.dump(obj, file)
    else:
        with open(path, 'wb') as file:
            pickle.dump(obj, file, protocol=pickle.HIGHEST_PROTOCOL)

def load(path):
    """ Loads a pickle, json or json.gz file.

    Parameters
    ----------
    path:
        The path to be loaded.

    Returns
    -------
    object
        The loaded object.
    """
    if path.endswith('.json.gz'):
        with gzip.open(path, 'r') as file:
            obj = json.loads(file.read().decode('utf-8'))
    elif path.endswith('.json'):
        with open(path,'r') as file:
            obj = json.load(file)
    else:
        with open(path, 'rb') as handle:
            obj = pickle.load(handle)
    return obj

def zip_file(path):
    """ Zips a file.

    Parameters
    ----------
    path:
        The path to the file.

    """
    with open(path, 'rb') as f_in:
        with gzip.open(path+'.gz', 'wb') as f_out:
            f_out.writelines(f_in)

def unzip_file(path):
    """ Unzips a .gz file.

    Parameters
    ----------
    path:
        The path to the .gz file.

    """
    assert path.endswith('.gz')
    with gzip.open(path, 'rb') as f_in:
        with open(path[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)



def download_url(url, out_path, log=True, chunk_size=10*1024*1024):
    """ Downloads a file from an url. If `out_path` is a directory, the file will be saved under the url basename.

    Parameters
    ----------
    url: str
        The url to be downloaded.
    out_path: str
        Path to save the downloaded file.
    log: bool, default True
        Whether to show a progress bar.
    chunk_size: int, default 10485760
        The chunk size of the download.

    """
    file_name = os.path.basename(url)
    if os.path.isdir(out_path) or out_path.endswith('/'):
        out_path += '/'+file_name
    r = requests.get(url, stream=True)
    r.raise_for_status()
    total = int(r.headers.get('content-length', 0))
    if log:
        print(f'Downloading {file_name}:')
        bar = tqdm(total=total, unit='iB', unit_scale=True, unit_divisor=chunk_size)
    with open(out_path, 'wb') as file:
        for data in r.iter_content(chunk_size=chunk_size):
            size = file.write(data)
            if log:
                bar.update(size)
    if log:
        bar.close()

def extract_tar(tar_path, out_path, extract_members=False):
    """ Extracts a tar file.

    Parameters
    ----------
    tar_path:
        The path to the tar file.
    out_path:
        The directory to extract to.
    extract_members: bool, default False
        If `True`, the tar file member will be directly extracted to `out_path`, instead of creating a subdirectory.
    """
    if extract_members:
        with tarfile.open(tar_path,'r') as file:
            for member in object.getmembers():
                file.extract(member, out_path)
    else:
        with tarfile.open(tar_path) as file:
            file.extractall(out_path)

def checkpoint(path):
    """ A decorator to checkpoint the result of a class method. It will check if the file specified by `path` exists, in which case the file is loaded and returned instead of running the decorated method. Otherwise the method is run and the result saved at the specified path.

    The path looks similar to a format string. The decorator will look in the attributes of the Owner of the decorated method to replace values in the format string. Example:

    .. code-block:: python

        @checkpoint('{root}/raw/{name}.pkl')
        def some_class_method(self):
            return {'hello': 'world'}

    This will replace `{root}` and `{name}` with `self.root` and `self.name` arguments of the class.

    """
    def decorator(function):
        def wrapper(self, *args, **kwargs):
            if os.path.exists(path.format(**self.__dict__)):
                return load(path.format(**self.__dict__))
            else:
                os.makedirs(os.path.dirname(path.format(**self.__dict__)), exist_ok=True)
                result = function(self, *args, **kwargs)
                save(result, path.format(**self.__dict__))
                return result
        return wrapper
    return decorator
