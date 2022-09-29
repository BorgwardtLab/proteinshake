import os
from tqdm import tqdm
import numpy as np

from proteinshake.utils import onehot

class Voxel():
    """ Voxel representation of a protein.

    Voxelizes a protein.

    Parameters
    ----------
    protein: dict
        A protein object.
    voxelsize: float
        The size of a voxel (in Angstrom).
    aggregation: str, defaul 'mean'
        How to aggregate labels of a voxel.

    """

    def __init__(self, protein, voxelsize, aggregation):
        resolution = 'atom' if 'atom' in protein else 'residue'
        self.protein = protein
        self.resolution = resolution
        labels = onehot(protein[resolution][f'{resolution}_type'])
        coords = np.stack([protein[resolution]['x'], protein[resolution]['y'], protein[resolution]['z']], axis=1)
        coords -= coords.min(axis=0) # translate to make all coords positive and flushed to the axes
        voxel_indices = (coords / voxelsize).astype(np.int32) # rasterize
        voxels = np.zeros(shape=np.concatenate([voxel_indices.max(axis=0)+1, labels.shape]))
        counts = np.zeros(shape=np.concatenate([voxel_indices.max(axis=0)+1, labels.shape]))
        voxel_indices = np.concatenate([voxel_indices, np.expand_dims(np.arange(len(coords)),1)], axis=1)
        voxels[tuple(voxel_indices.transpose())] = labels
        if aggregation == 'sum':
            self.voxel = voxels.sum(axis=-2)
        elif aggregation == 'mean':
            counts[tuple(voxel_indices.transpose())] = np.ones_like(labels)
            counts = counts.sum(axis=-2)
            self.voxel = np.divide(voxels.sum(axis=-2), counts, out=np.zeros_like(counts), where=counts!=0)



class VoxelDataset():
    """ Voxel representation of a protein structure dataset.

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
    voxelsize: float, default 10
        The size of a voxel (in Angstrom).
    aggregation: str, defaul 'mean'
        How to aggregate labels of a voxel.

    """

    def __init__(self, proteins, size, path, resolution='residue', voxelsize=10, aggregation='mean'):
        self.path = f'{path}/processed/voxel/{resolution}_voxelsize_{voxelsize}'
        self.voxels = (Voxel(protein, voxelsize, aggregation) for protein in proteins)
        self.size = size
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

    def torch(self, *args, **kwargs):
        from proteinshake.frameworks.torch import TorchVoxelDataset
        return TorchVoxelDataset(self.voxels, self.size, self.path+'.torch', *args, **kwargs)
