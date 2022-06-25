# -*- coding: utf-8 -*-
import os, gzip
import shutil
import os.path as osp
import pandas as pd
from tqdm import tqdm
import torch
from torch_geometric.utils import from_networkx
from torch_geometric.data import InMemoryDataset
from torch_geometric.data import extract_tar, download_url

from pdb_pyg.utils.pdb2graph import ProteinGraph

def pdb2pyg(*args, **kwargs):
    protein = ProteinGraph(*args, **kwargs)
    graph = protein.get_graph()
    graph = from_networkx(graph)
    graph.coord = graph.coord.float()
    if hasattr(graph, 'aa_idx'):
        graph.aa_idx = graph.aa_idx.float()
    graph.datapath = protein.datapath
    return graph
