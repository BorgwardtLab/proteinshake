import os
from itertools import repeat

import torch
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data, InMemoryDataset

class Dataset(InMemoryDataset):
    def __init__(self, root, data_list, transform=None, pre_transform=None, pre_filter=None):
        super().__init__(root, transform, pre_transform, pre_filter)
        if not os.path.exists(root):
            if self.pre_filter is not None:
                data_list = [data for data in data_list if self.pre_filter(data)]
            if self.pre_transform is not None:
                data_list = [self.pre_transform(data) for data in data_list]
            data, slices = self.collate(data_list)
            torch.save((data, slices), root)
        self.data, self.slices = torch.load(root)

    def get(self, idx):
        data = Data()
        for key in self.data.keys:
            item, slices = self.data[key], self.slices[key]
            if isinstance(item, torch.Tensor):
                s = list(repeat(slice(None), item.dim()))
                s[data.__cat_dim__(key, item)] = slice(slices[idx], slices[idx + 1])
                data[key] = item[s].clone()
            else:
                data[key] = [item[s] for s in range(slices[idx], slices[idx + 1])]
        return data
