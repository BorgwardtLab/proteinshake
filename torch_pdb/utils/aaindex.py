import os
import json
import numpy as np

aaindex_path = os.path.join(os.path.dirname(__file__), '..', 'pkg_data', 'aaindex.json')

PATH = os.path.realpath(os.path.dirname(__file__))

with open(aaindex_path, 'r') as file:
    aaindex = json.load(file)

def get_aaindex_features(amino_acid, format='array'):
    if format == 'array':
        return np.array(aaindex['properties'][amino_acid])
    elif format == 'dict':
        return {k:v for k,v in zip(aaindex['names'], aaindex['properties'][amino_acid])}
    else:
        raise 'Format not known.'
