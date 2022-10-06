import os
import itertools
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

    def __init__(self, protein, gridsize, voxelsize, aggregation):
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
            voxels = voxels.sum(axis=-2)
        elif aggregation == 'mean':
            counts[tuple(voxel_indices.transpose())] = np.ones_like(labels)
            counts = counts.sum(axis=-2)
            voxels = np.divide(voxels.sum(axis=-2), counts, out=np.zeros_like(counts), where=counts!=0)
        # pad to gridsize
        diff = gridsize-voxels.shape[:-1]
        paddings = np.stack([np.ceil(diff/2), np.floor(diff/2)],1).astype(np.int32)
        voxels = np.pad(voxels, (*paddings,(0,0)))
        self.voxel = voxels



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
    gridsize: tuple, default None
        The size of the grid in voxels as a 3-tuple of x,y,z edge lengths. If None (default), the dimensions of the largest protein in the dataset is used.
    voxelsize: float, default 10
        The size of a voxel (in Angstrom).
    aggregation: str, defaul 'mean'
        How to aggregate labels of a voxel.

    """

    def __init__(self, proteins, size, path, resolution='residue', gridsize=None, voxelsize=10, aggregation='mean'):
        self.path = f'{path}/processed/voxel/{resolution}_voxelsize_{voxelsize}'
        if gridsize is None:
            proteins, proteins_copy = itertools.tee(proteins)
            gridsize = np.array([[
                np.ptp(protein[resolution]['x']),
                np.ptp(protein[resolution]['y']),
                np.ptp(protein[resolution]['z'])
                ] for protein in proteins_copy]).max(0)
            gridsize = np.ceil(gridsize/voxelsize).astype(int)
        gridsize = np.array(gridsize)
        self.voxels = (Voxel(protein, gridsize, voxelsize, aggregation) for protein in proteins)
        self.size = size
        self.gridsize = gridsize
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

    def torch(self, *args, **kwargs):
        from proteinshake.frameworks.torch import TorchVoxelDataset
        return TorchVoxelDataset(self.voxels, self.size, self.path+'.torch', *args, **kwargs)
