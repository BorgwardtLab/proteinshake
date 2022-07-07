import os
import torch
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data
from tqdm import tqdm

from torch_pdb.utils import checkpoint, one_hot


class GraphDataset():

    def __init__(self, root, proteins, node_embedding=one_hot, eps=None, k=None, weighted_edges=False):
        assert not (eps is None and k is None), 'You must specify eps or k in the graph construction.'
        self.construction = 'k' if not k is None else 'eps'
        self.root = root
        self.k = k
        self.eps = eps
        self.weighted_edges = weighted_edges
        self.node_embedding = node_embedding
        self.name = f'emb_{node_embedding.__name__}'
        self.name += '_k_{k}' if self.construction == 'k' else f'_eps_{eps}'
        if self.weighted_edges:
            self.name += 'weighted'
        self.info = proteins
        self.proteins = self.convert(proteins)

    def protein2graph(self, protein):
        nodes = self.node_embedding(protein['sequence'])
        if self.construction == 'eps':
            mode = 'distance' if self.weighted_edges else 'connectivity'
            adj = radius_neighbors_graph(protein['coords'], radius=self.eps, mode=mode)
        elif self.construction == 'knn':
            adj = kneighbors_graph(protein['coords'], k=self.k)
        return (nodes, adj)

    def graph2pyg(self, graph, info={}):
        nodes = torch.Tensor(graph[0]).float()
        edges = from_scipy_sparse_matrix(graph[1])
        return Data(x=nodes, edge_index=edges[0].long(), edge_attr=edges[1].unsqueeze(1).float(), **info)

    @checkpoint('{root}/processed/graph/{name}.pkl')
    def convert(self, proteins):
        return [self.protein2graph(p) for p in tqdm(proteins, desc='Converting proteins to graphs')]

    @checkpoint('{root}/processed/graph/pyg__{name}.pkl')
    def pyg(self):
        return [self.graph2pyg(p, info=info) for p,info in zip(self.proteins,self.info)]

    @checkpoint('{root}/processed/graph/dgl__{name}.pkl')
    def dgl(self):
        raise NotImplementedError

    @checkpoint('{root}/processed/graph/nx__{name}.pkl')
    def nx(self):
        raise NotImplementedError
