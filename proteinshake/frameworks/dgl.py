import os
import dgl
import torch
from dgl.data import DGLDataset
from dgl import save_graphs, load_graphs
from tqdm import tqdm


class DGLGraphDataset(DGLDataset):
    """ Dataset class for graphs in DGL.

    Parameters
    ----------
    data_list: generator
        A generator of objects from a representation.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    transform: function
        A transform function to be applied in the __getitem__ method. Signature: transform(data, protein_dict) -> (data, protein_dict)
    """

    def __init__(self, data_list, size, path, transform=None):
        os.makedirs(path, exist_ok=True)
        self.path = path
        self.size = size
        self.transform = transform
        if not os.path.exists(f'{path}/{size-1}.pt'):
            for i, data_item in enumerate(tqdm(data_list, desc='Converting', total=size)):
                nodes, adj = data_item.data
                data = dgl.from_scipy(adj, eweight_name='edge_weight')
                if data_item.weighted_edges:
                    data.ndata[f'{data_item.resolution}'] = torch.tensor(nodes).long()
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
