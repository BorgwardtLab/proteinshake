import os

import dgl
import torch
from dgl.data import DGLDataset
from dgl import save_graphs, load_graphs

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


def graph2dgl(graph, info={}):
    g = dgl.from_scipy(graph[1])
    for key,value in info.items():
        if type(value) == list and len(value) == len(info['sequence']):
            try:
                g.ndata[key] = torch.tensor(value)
            except:
                pass
    return g

