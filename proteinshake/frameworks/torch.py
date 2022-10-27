import torch
import os
import numpy as np
from tqdm import tqdm
from torch.utils.data import Dataset


class TorchVoxelDataset(Dataset):
    """ Dataset class for voxels in torch.

    Parameters
    ----------
    data_list: generator
        A generator of objects from a representation.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    """

    def __init__(self, data_list, size, path, transform=None):
        os.makedirs(path, exist_ok=True)
        self.path = path
        self.size = size
        self.transform = transform
        if not os.path.exists(f'{path}/{size-1}.pt'):
            for i, data_item in enumerate(tqdm(data_list, desc='Converting', total=size)):
                data = torch.tensor(data_item.data).float().to_sparse()
                protein_dict = data_item.protein_dict
                torch.save((data, protein_dict), f'{path}/{i}.pt')

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        if idx > self.size - 1:
            raise StopIteration
        data, protein_dict = torch.load(f'{self.path}/{idx}.pt')
        data = data.to_dense()
        if not self.transform is None:
            data, protein_dict = self.transform(data, protein_dict)
        return data, protein_dict


class TorchPointDataset(Dataset):
    """ Dataset class for points in torch.

    Parameters
    ----------
    data_list: generator
        A generator of objects from a representation.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    """

    def __init__(self, data_list, size, path, transform=None):
        os.makedirs(path, exist_ok=True)
        self.path = path
        self.size = size
        self.transform = transform
        if not os.path.exists(f'{path}/{size-1}.pt'):
            for i, data_item in enumerate(tqdm(data_list, desc='Converting', total=size)):
                data = torch.tensor(data_item.data).float()
                protein_dict = data_item.protein_dict
                torch.save((data, protein_dict), f'{path}/{i}.pt')

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        if idx > self.size - 1:
            raise StopIteration
        data, protein_dict = torch.load(f'{self.path}/{idx}.pt')
        if not self.transform is None:
            data, protein_dict = self.transform(data, protein_dict)
        return data, protein_dict
