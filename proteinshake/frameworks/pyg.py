import torch
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data, Dataset as PygDataset
from proteinshake.frameworks.dataset import FrameworkDataset


class PygGraphDataset(FrameworkDataset, PygDataset):

    def convert_to_framework(self, data_item):
        nodes, adj = data_item.data
        return Data(
            x = torch.from_numpy(nodes),
            edge_index = from_scipy_sparse_matrix(adj)[0].long(),
            edge_attr = from_scipy_sparse_matrix(adj)[1].unsqueeze(1).float()
        )
