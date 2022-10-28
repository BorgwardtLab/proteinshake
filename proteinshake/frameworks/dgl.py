import dgl
import torch
from dgl.data import DGLDataset
from proteinshake.frameworks.dataset import FrameworkDataset


class DGLGraphDataset(FrameworkDataset, DGLDataset):

    def convert_to_framework(self, data_item):
        nodes, adj = data_item.data
        data = dgl.from_scipy(adj, eweight_name='edge_weight')
        if data_item.weighted_edges:
            data.ndata[f'{data_item.resolution}'] = torch.tensor(nodes).long()
        return data
