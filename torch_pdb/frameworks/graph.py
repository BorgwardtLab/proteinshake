import os
import torch
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data
from tqdm import tqdm

from torch_pdb.utils.io import load_if_exists


class GraphDataset():

    def __init__(self, root, proteins, node_embedding, eps=None, k=None, weighted_edges=False):
        assert not (eps is None and k is None), 'You must specify eps or k in the graph construction.'
        self.construction = 'k' if not k is None else 'eps'
        self.root = root
        self.k = k
        self.eps = eps
        self.weighted_edges = weighted_edges
        self.node_embedding = node_embedding
        self.name = f'k_{k}' if self.construction == 'k' else f'eps_{eps}'
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

    @load_if_exists('{root}/processed/{name}.graph')
    def convert(self, proteins):
        return [self.protein2graph(p) for p in tqdm(proteins, desc='Converting proteins to graphs')]

    @load_if_exists('{root}/processed/{name}.pyg')
    def pyg(self):
        return [self.graph2pyg(p, info=info) for p,info in zip(self.proteins,self.info)]

    @load_if_exists('{root}/processed/{name}.dgl')
    def dgl(self):
        raise NotImplementedError

    @load_if_exists('{root}/processed/{name}.nx')
    def nx(self):
        raise NotImplementedError
