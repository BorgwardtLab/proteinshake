import os
from tqdm import tqdm
import numpy as np

from proteinshake.utils import tokenize

class Point():
    """ Point representation of a protein.

    Parameters
    ----------
    protein: dict
        A protein object.

    """

    def __init__(self, protein):
        resolution = 'atom' if 'atom' in protein else 'residue'
        self.protein_dict = protein
        self.resolution = resolution
        labels = tokenize(protein[resolution][f'{resolution}_type'], resolution=resolution)
        coords =  np.stack([protein[resolution]['x'], protein[resolution]['y'], protein[resolution]['z']], axis=1)
        self.data = np.concatenate((coords,np.expand_dims(labels,-1)), axis=-1)



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

    def __init__(self, proteins, root, name, resolution='residue'):
        self.path = f'{root}/processed/point/{name}_{resolution}'
        self.points = (Point(protein) for protein in proteins)
        self.size = len(proteins)

    def torch(self, *args, **kwargs):
        from proteinshake.frameworks.torch import TorchPointDataset
        return TorchPointDataset(self.points, self.size, self.path+'.torch', *args, **kwargs)

    def tf(self, *args, **kwargs):
        from proteinshake.frameworks.tf import TensorflowPointDataset
        return TensorflowPointDataset(self.points, self.size, self.path+'.tf', *args, **kwargs)

    def np(self, *args, **kwargs):
        from proteinshake.frameworks.np import NumpyPointDataset
        return NumpyPointDataset(self.points, self.size, self.path+".np", *args, **kwargs)
