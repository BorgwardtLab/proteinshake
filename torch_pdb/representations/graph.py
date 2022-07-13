import os
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from tqdm import tqdm
import numpy as np

from torch_pdb.utils import checkpoint, one_hot


class GraphDataset():
    """ Graph representation of a protein structure dataset.

    Converts a protein object to a graph by using a k-nearest-neighbor or epsilon-neighborhood approach. Define either `k` or `eps` to determine which one is used.

    Also embeds the protein sequence to attributes of the graph nodes using a supplied embedding function. See `utils.embeddings` for examples. If a list of functions is passed to `embedding`, the resulting features will be concatenated.

    Parameters
    ----------
    embedding: Union[function, list]
        A function or list of functions for embedding the protein sequence to node attributes.
    eps: float
        The epsilon radius to be used in graph construction (in Angstrom).
    k: int
        The number of neighbors to be used in the k-NN graph.
    weighted_edges: bool, default False
        If `True`, edges are attributed with their euclidean distance. If `False`, edges are unweighted.

    """

    def __init__(self, root, proteins, embedding=one_hot, eps=None, k=None, weighted_edges=False):
        assert not (eps is None and k is None), 'You must specify eps or k in the graph construction.'
        self.construction = 'knn' if not k is None else 'eps'
        self.root = root
        self.k = k
        self.eps = eps
        self.weighted_edges = weighted_edges
        self.embedding = embedding
        self.name = f'emb_{embedding.__name__}'
        self.name += '_k_{k}' if self.construction == 'knn' else f'_eps_{eps}'
        if self.weighted_edges:
            self.name += 'weighted'
        self.info = proteins
        self.proteins = self.convert(proteins)

    def protein2graph(self, protein):
        nodes = self.embedding(protein['sequence'])
        mode = 'distance' if self.weighted_edges else 'connectivity'
        if self.construction == 'eps':
            adj = radius_neighbors_graph(protein['coords'], radius=self.eps, mode=mode)
        elif self.construction == 'knn':
            adj = kneighbors_graph(protein['coords'], n_neighbors=self.k, mode=mode)
        return (nodes, adj)

    @checkpoint('{root}/processed/graph/{name}.pkl')
    def convert(self, proteins):
        return [self.protein2graph(p) for p in tqdm(proteins, desc='Converting proteins to graphs')]

    @checkpoint('{root}/processed/graph/{name}.pyg.pkl')
    def pyg(self):
        import torch
        from torch_geometric.utils import from_scipy_sparse_matrix
        from torch_geometric.data import Data
        def graph2pyg(graph, info={}):
            nodes = torch.Tensor(graph[0]).float()
            edges = from_scipy_sparse_matrix(graph[1])
            return Data(x=nodes, edge_index=edges[0].long(), edge_attr=edges[1].unsqueeze(1).float(), **info)
        return [graph2pyg(p, info=info) for p,info in zip(self.proteins,self.info)]

    def dgl(self):
        from dgl.data import DGLDataset
        from dgl import save_graphs, load_graphs
        import dgl
        import torch
        def graph2dgl(graph, info={}):
            g = dgl.from_scipy(graph[1])
            for key,value in info.items():
                if type(value) == list and len(value) == len(info['sequence']):
                    try:
                        g.ndata[key] = torch.tensor(value)
                    except:
                        pass
            return g
        class Dataset(DGLDataset):
            def __init__(self, path, proteins, info):
                if os.path.exists(path):
                    self.proteins = load_graphs(path)
                else:
                    self.proteins = [graph2dgl(p, info=info) for p,info in zip(proteins,info)]
                    save_graphs(path, self.proteins)
            def __getitem__(self, i):
                return self.proteins[i]
            def __len__(self):
                return len(self.proteins)
        ds = Dataset(f'{self.root}/processed/graph/{self.name}.dgl.pkl',  self.proteins, self.info)
        return ds

    @checkpoint('{root}/processed/graph/{name}.nx.pkl')
    def nx(self):
        import networkx as nx
        def graph2nx(graph, info={}):
            g = nx.from_scipy_sparse_matrix(graph[1])
            g.add_nodes_from(graph[0])
            for key,value in info.items():
                if type(value) == list and len(value) == len(info['sequence']):
                    nx.set_node_attributes(g, value, key)
                else:
                    g.graph[key] = value
            return g
        return [graph2nx(p, info=info) for p,info in zip(self.proteins,self.info)]
