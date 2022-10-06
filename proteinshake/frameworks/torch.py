import torch
import os
import numpy as np
from tqdm import tqdm
from torch.utils.data import Dataset as TorchDataset

def convert_to_tensor(value):
    try:
        return torch.tensor(value).float()
    except:
        return value

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

    def __init__(self, voxels, size, path, transform=None):
        self.size = size
        self.transform = transform
        if not os.path.exists(path):
            _voxels, _protein_dicts = [],[]
            for voxel in tqdm(voxels, desc='Converting to voxel', total=size):
                _voxels.append(torch.tensor(voxel.voxel).float().to_sparse())
                _protein_dicts.append({
                    'protein': {k:convert_to_tensor(v) for k,v in voxel.protein['protein'].items()},
                    voxel.resolution: {k:convert_to_tensor(v) for k,v in voxel.protein[voxel.resolution].items() if k not in ['atom_type','atom_number','residue_type','residue_number','x','y','z']},
                })
            torch.save((_voxels, _protein_dicts), path)
            del voxels
        self.voxels, self.protein_dicts = torch.load(path)

    def __len__(self):
        return len(self.voxels)

    def __getitem__(self, idx):
        if not self.transform is None:
            return self.transform(self.voxels[idx].to_dense(), self.protein_dicts[idx])
        return self.voxels[idx].to_dense(), self.protein_dicts[idx]



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

    def __init__(self, points, size, path, transform=None):
        self.size = size
        self.transform = transform
        if not os.path.exists(path):
            coords, labels, protein_dicts = [], [], []
            for point in tqdm(points, desc='Converting to point', total=size):
                coords.append(torch.tensor(point.coords))
                labels.append(torch.tensor(point.labels))
                protein_dicts.append({
                    'protein': {k:convert_to_tensor(v) for k,v in point.protein['protein'].items()},
                    point.resolution: {k:convert_to_tensor(v) for k,v in point.protein[point.resolution].items() if k not in ['atom_type','atom_number','residue_type','residue_number','x','y','z']},
                })
            torch.save((coords, labels, protein_dicts), path)
            del coords, labels
        self.coords, self.labels, self.protein_dicts = torch.load(path)

    def __len__(self):
        return len(self.coords)

    def __getitem__(self, idx):
        if not self.transform is None:
            return self.transform(self.coords[idx], self.labels[idx], self.protein_dicts[idx])
        return self.coords[idx], self.labels[idx], self.protein_dicts[idx]
