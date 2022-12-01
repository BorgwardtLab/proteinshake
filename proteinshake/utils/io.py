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
import re
from joblib import Parallel

import pandas as pd
from tqdm import tqdm
from fastavro import writer as avro_writer, reader as avro_reader, parse_schema as parse_avro_schema

AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
AA_ONE_TO_THREE = {v:k for k, v in AA_THREE_TO_ONE.items()}

class ProgressParallel(Parallel):
    """ Extends joblib's Parallel with a progress bar.

    Parameters
    ----------
    total:
        The total number of jobs.
    desc: A description to display.
    """
    def __init__(self, total=None, desc=None, *args, **kwargs):
        self._total = total
        self._desc = desc
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        with tqdm(total=self._total, desc=self._desc) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()

def fx2str(fx):
    """ Converts a function to a string representation.

    Parameters
    ----------
    fx: function
        A function.

    Returns
    -------
    str
        The stringified function.
    """
    return re.sub('(<.*?)\\s.*(>)', r'\1\2', fx.__repr__())

def avro_schema_from_protein(protein):
    """ Guesses the avro schema from a dictionary.

    Parameters
    ----------
    protein: dict
        A protein dictionary.

    Returns
    -------
    schema
        An avro schema.
    """
    typedict = {'int':'int', 'float':'float', 'str':'string', 'bool':'boolean'}
    def field_spec(k,v):
        if type(v) == dict:
            return {'name':k, 'type':{'name':k, 'type':'record', 'fields': [field_spec(_k,_v) for _k,_v in v.items()]}}
        elif type(v) == list:
            return {'name':k, 'type':{'type': 'array', 'items': typedict[type(v[0]).__name__]}}
        elif type(v).__name__ in typedict:
            return {'name':k, 'type': typedict[type(v).__name__]}
        else:
            raise f'All fields in a protein object need to be either int, float, bool or string, not {type(v).__name__}'
    schema = {
        'name': 'Protein',
        'namespace': 'Dataset',
        'type': 'record',
        'fields': [field_spec(k,v) for k,v in protein.items()],
    }
    return parse_avro_schema(schema)

def write_avro(proteins, path):
    """ Writes a list of protein dictionaries to an avro file.

    Parameters
    ----------
    proteins: list
        The list of proteins.
    path:
        The path to the output file.
    """
    schema = avro_schema_from_protein(proteins[0])
    with open(path, 'wb') as file:
        avro_writer(file, schema, proteins, metadata={'number_of_proteins':str(len(proteins))})

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
    os.remove(path)



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
            for member in tqdm(file.getmembers(), desc='Extracting', total=len(file.getmembers())):
                file.extract(member, out_path)
    else:
        with tarfile.open(tar_path) as file:
            file.extractall(out_path, members=tqdm(file, desc='Extracting', total=len(file.getmembers())))

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

def protein_to_pdb(protein, path):
    """ Write coordinate list from atom dict to a PDB file.

    Parameters
    ------------
    protein:
        protein data dictionary. must be at the 'atom' resolution.
    path:
        Path to write PDB file.


    """
    # ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
    df = pd.DataFrame(protein['atom'])
    df['residue_name_full'] = df['residue_type'].apply(lambda x : AA_ONE_TO_THREE[x])
    if 'chain_id' not in df.columns:
        df['chain_id'] = ['A'] * len(df)
    df['occupancy'] = [1.00] * len(df)
    df['temp'] = [20.00] * len(df)
    df['element'] = df['atom_type'].apply(lambda x: x[:1])
    lines = []
    for row in df.itertuples():
        line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   " \
               "{:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          " \
               "{:>2s}{:2s}".format(
                                    'ATOM',
                                    row.atom_number,
                                    row.atom_type,
                                    ' ',
                                    row.residue_name_full,
                                    row.chain_id,
                                    row.residue_number,
                                    ' ',
                                    row.x,
                                    row.y,
                                    row.z,
                                    row.occupancy,
                                    row.temp,
                                    row.element,
                                    '  '
                                 )
        lines.append(line)
    with open(path, "w") as p:
        p.write("\n".join(lines))
    pass
