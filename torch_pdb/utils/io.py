import os
import tarfile
import torch
import pickle
import h5py
import requests
from tqdm import tqdm

# decorator to return saved file if it exists
def checkpoint(path):
    def decorator(function):
        def wrapper(self, *args, **kwargs):
            if os.path.exists(path.format(**self.__dict__)):
                return torch.load(path.format(**self.__dict__))
            else:
                os.makedirs(os.path.dirname(path.format(**self.__dict__)), exist_ok=True)
                result = function(self, *args, **kwargs)
                torch.save(result, path.format(**self.__dict__))
                return result
        return wrapper
    return decorator

def save(obj, path):
    if path.endswith('.h5'):
        with h5py.File(path, 'w') as file:
            file.create_dataset('data', data=obj)
    else:
        with open(path, 'wb') as file:
            pickle.dump(obj, file, protocol=pickle.HIGHEST_PROTOCOL)

def load(path):
    if path.endswith('.h5'):
        with h5py.File(path, 'r') as file:
            obj = file['data']
    else:
        with open(path, 'rb') as handle:
            obj = pickle.load(handle)
    return obj

def download_url(url, out_path, log=True, chunk_size=10*1024):
    file_name = os.path.basename(url)
    if os.path.isdir(out_path):
        out_path += '/'+file_name
    r = requests.get(url, stream=True)
    total = int(r.headers.get('content-length', 0))
    if log:
        print(f'Downloading {file_name}:')
        bar = tqdm(total=total, unit='iB', unit_scale=True, unit_divisor=chunk_size)
    with open(out_path, 'wb') as file:
        for data in r.iter_content(chunk_size=chunk_size):
            size = file.write(data)
            if log:
                bar.update(size)

def extract_tar(tar_path, out_path, extract_members=False):
    if replace:
        with tarfile.open(tar_path,'r') as file:
            for member in object.getmembers():
                file.extract(member, out_path)
    else:
        with tarfile.open(tar_path) as file:
            file.extractall(out_path)
