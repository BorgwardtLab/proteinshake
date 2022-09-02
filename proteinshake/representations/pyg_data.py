import os
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
