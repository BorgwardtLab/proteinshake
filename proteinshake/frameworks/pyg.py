import os
from itertools import repeat
import torch
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data, InMemoryDataset
from proteinshake.utils import fx2str
from tqdm import tqdm


def info2pyg(info):
    """ Try to convert info dictionary to tensor."""
    new_info = {}
    for k,v in info.items():
        if k in ['x','y','z']:
            continue
        try:
            new_info[k] = torch.Tensor(v).long()
        except:
            new_info[k] = v
            pass
    return new_info

class PygGraphDataset(InMemoryDataset):
    """ Dataset class for graph in pytorch-geometric.

    Parameters
    ----------
    graphs: generator
        A generator of graph objects from GraphDataset.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    transform: function
        Compare torch_geometric.data.Dataset.
    pre_transform: function
        Compare torch_geometric.data.Dataset.
    pre_filter: function
        Compare torch_geometric.data.Dataset.
    """
    def __init__(self, graphs, size, path, transform=None, pre_transform=None, pre_filter=None):
        super().__init__(None, transform, pre_transform, pre_filter)
        self.size = size
        transforms_repr = fx2str(transform)+fx2str(pre_transform)+fx2str(pre_filter)
        if not os.path.exists(path):
            data_list = [Data(
                x = torch.from_numpy(graph.nodes),
                edge_index = from_scipy_sparse_matrix(graph.adj)[0].long(),
                edge_attr = from_scipy_sparse_matrix(graph.adj)[1].unsqueeze(1).float(),
                **info2pyg(graph.protein['protein']),
                **info2pyg(graph.protein[graph.resolution])
            ) for graph in tqdm(graphs, desc='Converting to graph', total=size)]
            if self.pre_filter is not None:
                data_list = [data for data in data_list if self.pre_filter(data)]
            if self.pre_transform is not None:
                data_list = [self.pre_transform(data) for data in data_list]
            data, slices = self.collate(data_list)
            torch.save((transforms_repr, data, slices), path)
            del data_list, data, slices
        original_repr, self.data, self.slices = torch.load(path)
        assert original_repr == transforms_repr, f'The transforms are not the same as when the dataset was created. If you want to change them, delete the file at {path}'

    def __len__(self):
        return self.size

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
