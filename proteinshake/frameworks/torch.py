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
            voxels = [voxel.voxel for voxel in tqdm(voxels, desc='Converting to voxel', total=size)]
            max_dims = np.max([v.shape[:-1] for v in voxels], axis=0)
            def pad(voxel):
                diff = max_dims-voxel.shape[:-1]
                paddings = np.stack([np.ceil(diff/2), np.floor(diff/2)],1).astype(np.int32)
                padded = np.pad(voxel, (*paddings,(0,0)))
                return padded
            voxels = np.stack([pad(v) for v in voxels])
            voxels = torch.tensor(voxels).to_sparse()
            torch.save(voxels, path)
            del voxels
        self.voxels = torch.load(path)


    def __len__(self):
        return len(self.voxels)

    def __getitem__(self, idx):
        return self.voxels[idx].to_dense()
