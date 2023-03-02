import numpy as np
from scipy.spatial.distance import pdist, squareform

def global_distance_test(coordsA, coordsB):
    """
    Computes GDT. Input: aligned coordinates of both structures.
    """
    dist = np.sqrt(((coordsA-coordsB)**2).sum(1))
    thresholds = [1,2,4,8]
    return np.mean([np.mean(dist<=t) for t in thresholds])

def local_distance_difference_test(coordsA, coordsB):
    # parameters from 
    R0 = 15
    thresholds = [0.5,1,2,4]

    distA = squareform(pdist(coordsA))
    np.fill_diagonal(distA, np.nan)
    distA[distA>R0] = np.nan

    distB = squareform(pdist(coordsB))
    np.fill_diagonal(distB, np.nan)
    distB[distB>R0] = np.nan

    diff = np.abs(distA.reshape(-1) - distB.reshape(-1))
    return np.mean([np.sum(diff<=t)/np.isfinite(diff).sum() for t in thresholds])


