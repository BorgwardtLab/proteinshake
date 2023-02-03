from scipy.spatial.distance import cdist
import numpy as np

def gdt(superposition):
    df = superposition.df['ATOM']
    df = df[df['atom_name'] == 'CA']
    A = df[df['chain_id'] == 'A']
    B = df[df['chain_id'] == 'B']
    coordsA = np.array(list(zip(A['x_coord'], A['y_coord'], A['z_coord'])))
    coordsB = np.array(list(zip(B['x_coord'], B['y_coord'], B['z_coord'])))
    dist = np.sqrt(((coordsA-coordsB)**2).sum(1))
    thresholds = [1,2,4,8]
    return np.mean([np.mean(dist<=t) for t in thresholds])
