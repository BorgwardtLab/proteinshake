import os
import json
import numpy as np

PATH = os.path.realpath(os.path.dirname(__file__))

with open('{}/_aaindex/aaindex.json'.format(PATH),'r') as file:
    aaindex = json.load(file)

def get_aaindex_features(amino_acid, format='array'):
    if format == 'array':
        return np.array(aaindex['properties'][amino_acid])
    elif format == 'dict':
        return {k:v for k,v in zip(aaindex['names'], aaindex['properties'][amino_acid])}
    else:
        raise 'Format not known.'
