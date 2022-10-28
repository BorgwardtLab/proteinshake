import networkx as nx
from proteinshake.frameworks.dataset import FrameworkDataset


class NetworkxGraphDataset(FrameworkDataset):

    def convert_to_framework(self, data_item):
        nodes, adj = data_item.data
        data = nx.from_scipy_sparse_array(adj)
        data.add_nodes_from(nodes)
        return data
