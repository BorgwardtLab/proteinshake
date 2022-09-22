import os
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from tqdm import tqdm
import numpy as np

from proteinshake.utils import tokenize


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

    def __init__(self, dataset, resolution='residue', eps=None, k=None, weighted_edges=False):
        assert not (eps is None and k is None), 'You must specify eps or k in the graph construction.'
        self.construction = 'knn' if not k is None else 'eps'
        self.resolution = resolution
        self.weighted_edges = weighted_edges
        self.eps = eps
        self.k = k
        self.name = f'resolution_{resolution}_'
        self.name += f'knn_{k}' if self.construction == 'knn' else f'eps_{eps}'
        if self.weighted_edges:
            self.name += '_weighted'
        self.dataset = dataset
        os.makedirs(f'{dataset.root}/processed/graph', exist_ok=True)
        dataset.download_precomputed(resolution=resolution)

    def convert(self, protein):
        resolution = 'atom' if 'atom' in protein else 'residue'
        nodes = tokenize(protein['protein']['sequence'])
        mode = 'distance' if self.weighted_edges else 'connectivity'
        coords = np.stack([protein[resolution]['x'], protein[resolution]['y'], protein[resolution]['z']], axis=1)
        if self.construction == 'eps':
            adj = radius_neighbors_graph(coords, radius=self.eps, mode=mode)
        elif self.construction == 'knn':
            n_neighbors = min(len(coords) - 1, self.k) # reduce k if protein is smaller than self.k
            adj = kneighbors_graph(coords,  n_neighbors=n_neighbors, mode=mode)
        return nodes, adj

    def pyg(self, *args, **kwargs):
        from proteinshake.frameworks import PygGraphDataset
        return PygGraphDataset(self, *args, **kwargs)
