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
import pandas as pd
import numpy as np
from pathlib import Path
from tqdm import tqdm
from fastavro import writer as avro_writer, reader as avro_reader, parse_schema as parse_avro_schema

AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
AA_ONE_TO_THREE = {v:k for k, v in AA_THREE_TO_ONE.items()}

class Generator(object):
    def __init__(self, generator, length):
        self.generator = generator
        self.length = length

    def __len__(self): 
        return self.length

    def __iter__(self):
        return self.generator

    def __next__(self):
        return next(self.generator)

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
            return {'name':k, 'type':{'type': 'array', 'items': typedict[type(v[0]).__name__] if len(v)>0 else 'string'}}
        elif type(v).__name__ in typedict:
            return {'name':k, 'type': typedict[type(v).__name__]}
        else:
            raise TypeError(f"All fields in a protein object need to be either int, float, bool or string, not {type(v).__name__}")
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
    path = Path(path)
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
        with gzip.open(Path(path), 'w') as file:
            file.write(json.dumps(obj).encode('utf-8'))
    elif path.endswith('.json'):
        with open(Path(path),'w') as file:
            json.dump(obj, file)
    elif path.endswith('.npy'):
        np.save(path, obj)
    else:
        with open(Path(path), 'wb') as file:
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
        with gzip.open(Path(path), 'r') as file:
            obj = json.loads(file.read().decode('utf-8'))
    elif path.endswith('.json'):
        with open(Path(path),'r') as file:
            obj = json.load(file)
    elif path.endswith('.npy'):
        obj = np.load(path)
    else:
        with open(Path(path), 'rb') as handle:
            obj = pickle.load(handle)
    return obj

def zip_file(path):
    """ Zips a file.

    Parameters
    ----------
    path:
        The path to the file.

    """
    with open(Path(path), 'rb') as f_in:
        with gzip.open(Path(path+'.gz'), 'wb') as f_out:
            f_out.writelines(f_in)
    return path+'.gz'

def unzip_file(path, remove=True):
    """ Unzips a .gz file.

    Parameters
    ----------
    path:
        The path to the .gz file.

    """
    assert path.endswith('.gz')
    with gzip.open(Path(path), 'rb') as f_in:
        with open(Path(path[:-3]), 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    if remove:
        os.remove(path)
    return path[:-3]



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
    out_path = Path(out_path)
    with open(out_path, 'wb') as file:
        for data in r.iter_content(chunk_size=chunk_size):
            size = file.write(data)
            if log:
                bar.update(size)
    if log:
        bar.close()

def extract_tar(tar_path, out_path, extract_members=False, strip=0):
    """ Extracts a tar file.

    Parameters
    ----------
    tar_path:
        The path to the tar file.
    out_path:
        The directory to extract to.
    extract_members: bool, default False
        If `True`, the tar file member will be directly extracted to `out_path`, instead of creating a subdirectory.
    strip: int, default 0
        Remove `strip` folder hierarchies from the path of the extracted file.
    """
    def get_members(file):
        for member in file.getmembers():
            parts = Path(member.path).parts
            member.path = Path(*parts[min(strip, len(parts)-1):])
            yield member

    out_path = Path(out_path)
    with tarfile.open(tar_path,'r') as file:
        len_members = len(file.getmembers())
        members = get_members(file)
        if extract_members:
            for member in tqdm(members, desc='Extracting', total=len_members):
                file.extract(member, out_path)
        else:
            file.extractall(out_path, members=tqdm(file, desc='Extracting', total=len_members))

def protein_to_pdb(protein, path):
    """ Write coordinate list from atom dict to a PDB file.

    Parameters
    ------------
    protein:
        protein data dictionary. must be at the 'atom' resolution.
    path:
        Path to write PDB file.


    """
    path = Path(path)
    # ATOM      1  N   PRO A   1       8.316  21.206  21.530  1.00 17.44           N
    try:
        df = pd.DataFrame(protein['atom'])
        mode = 'atom'
    except KeyError:
        df = pd.DataFrame(protein['residue'])
        mode = 'residue'

    df['residue_name_full'] = df['residue_type'].apply(lambda x : AA_ONE_TO_THREE[x])
    if 'chain_id' not in df.columns:
        df['chain_id'] = 'A'
    df['occupancy'] = [1.00] * len(df)
    df['temp'] = [20.00] * len(df)

    if mode == 'atom':
        df['element'] = df['atom_type'].apply(lambda x: x[:1])
    else:
        df['element'] = 'C'
        df['atom_number'] = df['residue_number']
        df['atom_type'] = 'CA'

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
