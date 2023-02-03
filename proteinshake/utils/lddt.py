from scipy.spatial.distance import pdist
import numpy as np

def lddt(proteinA, proteinB, resolution='residue'):
    coordsA = list(zip(proteinA[resolution]['x'],proteinA[resolution]['y'],proteinA[resolution]['z']))
    coordsB = list(zip(proteinB[resolution]['x'],proteinB[resolution]['y'],proteinB[resolution]['z']))
    dA = pdist(coordsA)[np.triu_indices(len(coordsA))]
    dB = pdist(coordsB)[np.triu_indices(len(coordsB))]
    diff = np.absolute(dA-dB)
    thresholds = [0.5,1,2,4] # from Mariani et al 2013
    return np.mean([np.mean(diff<=t) for t in thresholds])
