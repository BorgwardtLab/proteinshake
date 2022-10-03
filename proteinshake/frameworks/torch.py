import torch
import os
import numpy as np
from tqdm import tqdm
from torch.utils.data import Dataset as TorchDataset

class TorchVoxelDataset(TorchDataset):
    """ Dataset class for voxel in torch.

    Parameters
    ----------
    voxels: generator
        A generator of voxel objects from VoxelDataset.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    """

    def __init__(self, voxels, size, path):
        self.size = size
        if not os.path.exists(path):
            voxels = [torch.tensor(voxel.voxel).to_sparse() for voxel in tqdm(voxels, desc='Converting to voxel', total=size)]
            torch.save(voxels, path)
            del voxels
        self.voxels = torch.load(path)


    def __len__(self):
        return len(self.voxels)

    def __getitem__(self, idx):
        return self.voxels[idx].to_dense()


class TorchPointDataset(TorchDataset):
    """ Dataset class for point in torch.

    Parameters
    ----------
    points: generator
        A generator of point objects from PointDataset.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    """

    def __init__(self, points, size, path):
        self.size = size
        if not os.path.exists(path):
            points = [torch.tensor(point.coords) for point in tqdm(points, desc='Converting to point', total=size)]
            torch.save(points, path)
            del points
        self.points = torch.load(path)


    def __len__(self):
        return len(self.points)

    def __getitem__(self, idx):
        return self.points[idx]
