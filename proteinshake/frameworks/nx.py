import os
import networkx as nx
from tqdm import tqdm
from proteinshake.utils import save, load


class NetworkxGraphDataset():
    """ Dataset class for graphs in NetworkX.

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
        if not os.path.exists(f'{path}/{size-1}.pkl'):
            for i, data_item in enumerate(tqdm(data_list, desc='Converting', total=size)):
                nodes, adj = data_item.data
                data = nx.from_scipy_sparse_matrix(adj)
                data.add_nodes_from(nodes)
                protein_dict = data_item.protein_dict
                save((data, protein_dict), f'{path}/{i}.pkl')

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        if idx > self.size - 1:
            raise StopIteration
        data, protein_dict = load(f'{self.path}/{idx}.pkl')
        if not self.transform is None:
            data, protein_dict = self.transform(data, protein_dict)
        return data, protein_dict
