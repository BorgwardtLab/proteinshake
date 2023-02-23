import os
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from tqdm import tqdm
import numpy as np

from proteinshake.utils import tokenize

class Graph():
    """ Graph representation of a protein.

    Converts a protein object to a graph by using a k-nearest-neighbor or epsilon-neighborhood approach. Define either `k` or `eps` to determine which one is used.

    Parameters
    ----------
    protein: dict
        A protein object.
    construction: str
        Whether to use knn or eps construction.
    eps: float
        The epsilon radius to be used in graph construction (in Angstrom).
    k: int
        The number of neighbors to be used in the k-NN graph.
    weighted_edges: bool, default False
        If `True`, edges are attributed with their euclidean distance. If `False`, edges are unweighted.

    """

    def __init__(self, protein, construction, k, eps, weighted_edges):
        resolution = 'atom' if 'atom' in protein else 'residue'
        mode = 'distance' if weighted_edges else 'connectivity'
        coords = np.stack([protein[resolution]['x'], protein[resolution]['y'], protein[resolution]['z']], axis=1)
        nodes = tokenize(protein[resolution][f'{resolution}_type'], resolution=resolution)
        if construction == 'eps':
            adj = radius_neighbors_graph(coords, radius=eps, mode=mode)
        elif construction == 'knn':
            n_neighbors = min(len(coords) - 1, k) # reduce k if protein is smaller than self.k
            adj = kneighbors_graph(coords,  n_neighbors=n_neighbors, mode=mode)
        self.protein_dict = protein
        self.resolution = resolution
        self.data = (nodes, adj)
        self.weighted_edges = weighted_edges



class GraphDataset():
    """ Graph representation of a protein structure dataset.

    Parameters
    ----------
    proteins: generator
        A generator of protein objects from a Dataset.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    resolution: str, default 'residue'
        Resolution of the proteins to use in the graph representation. Can be 'atom' or 'residue'.
    eps: float
        The epsilon radius to be used in graph construction (in Angstrom).
    k: int
        The number of neighbors to be used in the k-NN graph.
    weighted_edges: bool, default False
        If `True`, edges are attributed with their euclidean distance. If `False`, edges are unweighted.

    """

    def __init__(self, proteins, root, name, resolution='residue', eps=None, k=None, weighted_edges=False):
        assert not (eps is None and k is None), 'You must specify eps or k in the graph construction.'
        construction = 'knn' if not k is None else 'eps'
        param = k if construction == 'knn' else eps
        weighted = '_weighted' if weighted_edges else ''
        self.path = f'{root}/processed/graph/{name}_{resolution}_{construction}_{param}{weighted}'
        self.graphs = (Graph(protein, construction, k, eps, weighted_edges) for protein in proteins)
        self.size = len(proteins)
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

    def pyg(self, *args, **kwargs):
        from proteinshake.frameworks.pyg import PygGraphDataset
        return PygGraphDataset(self.graphs, self.size, self.path+'.pyg', *args, **kwargs)

    def dgl(self, *args, **kwargs):
        from proteinshake.frameworks.dgl import DGLGraphDataset
        return DGLGraphDataset(self.graphs, self.size, self.path+'.dgl', *args, **kwargs)

    def nx(self, *args, **kwargs):
        from proteinshake.frameworks.nx import NetworkxGraphDataset
        return NetworkxGraphDataset(self.graphs, self.size, self.path+'.nx', *args, **kwargs)
