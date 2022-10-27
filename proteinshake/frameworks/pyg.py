import os
from itertools import repeat
import torch
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data, Dataset
from proteinshake.utils import fx2str
from tqdm import tqdm


class PygGraphDataset(Dataset):
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
    def __init__(self, data_list, size, path, transform=None, pre_transform=None, pre_filter=None):
        super().__init__(None, transform, pre_transform, pre_filter)
        os.makedirs(path, exist_ok=True)
        self.size = size
        self.path = path
        self.transform = transform
        transforms_repr = fx2str(pre_transform)+fx2str(pre_filter)
        if not os.path.exists(f'{path}/{size-1}.pt'):
            torch.save(transforms_repr, f'{path}/transforms.pt')
            for i, data_item in enumerate(tqdm(data_list, desc='Converting', total=size)):
                protein_dict = data_item.protein_dict
                nodes, adj = data_item.data
                data = Data(
                    x = torch.from_numpy(nodes),
                    edge_index = from_scipy_sparse_matrix(adj)[0].long(),
                    edge_attr = from_scipy_sparse_matrix(adj)[1].unsqueeze(1).float()
                )
                if self.pre_filter is not None:
                    data, protein_dict = self.pre_filter(data, protein_dict)
                if self.pre_transform is not None:
                    data, protein_dict = self.pre_transform(data, protein_dict)
                torch.save((data, protein_dict), f'{path}/{i}.pt')
        original_repr = torch.load(f'{path}/transforms.pt')
        assert original_repr == transforms_repr, f'The transforms are not the same as when the dataset was created. If you want to change them, delete the folder at {path}'

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        if idx > self.size - 1:
            raise StopIteration
        data, protein_dict = torch.load(f'{self.path}/{idx}.pt')
        if not self.transform is None:
            return self.transform(data, protein_dict)
        else:
            return data, protein_dict
