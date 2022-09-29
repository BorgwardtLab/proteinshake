import os
from tqdm import tqdm
import numpy as np

from proteinshake.utils import onehot

class Point():
    """ Point representation of a protein.

    Parameters
    ----------
    protein: dict
        A protein object.

    """

    def __init__(self, protein):
        resolution = 'atom' if 'atom' in protein else 'residue'
        self.protein = protein
        self.resolution = resolution
        self.coords =  np.stack([protein[resolution]['x'], protein[resolution]['y'], protein[resolution]['z']], axis=1)



class PointDataset():
    """ Point representation of a protein structure dataset.

    Parameters
    ----------
    proteins: generator
        A generator of protein objects from a Dataset.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    resolution: str, default 'residue'
        Resolution of the proteins to use in the graph representation. Can be 'atom' or 'residue'.

    """

    def __init__(self, proteins, size, path, resolution='residue'):
        self.path = f'{path}/processed/point/{resolution}'
        self.points = (Point(protein) for protein in proteins)
        self.size = size
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

    def torch(self, *args, **kwargs):
        from proteinshake.frameworks.torch import TorchPointDataset
        return TorchPointDataset(self.points, self.size, self.path+'.torch', *args, **kwargs)
