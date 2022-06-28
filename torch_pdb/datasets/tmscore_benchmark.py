# -*- coding: utf-8 -*-
import os, gzip
import shutil
import os.path as osp

import wget
import pandas as pd
from tqdm import tqdm
import torch
from torch_geometric.data import download_url, extract_tar
from torch_geometric.utils import from_networkx
from torch_geometric.data import InMemoryDataset
from biopandas.pdb import PandasPdb

from torch_pdb.utils.convert import pdb2pyg

class TMScoreBenchmark(InMemoryDataset):
    def __init__(self, name, root='/tmp/var', transform=None, pre_transform=None,
                 url="https://github.com/BorgwardtLab/torch-pdb/releases/download/v1.0.0/tm-bench.tar.gz",
                 **kwargs):
        '''
            - name (str): name of the dataset
            - root (str): root directory to store the dataset folder
            - transform, pre_transform (optional): transform/pre-transform graph objects
            - meta_dict: dictionary that stores all the meta-information about data. Default is None,
                    but when something is passed, it uses its information. Useful for debugging for external contributers.
        '''
        # filename="PDBbind_v2020_refined.tar.gz",
        self.name = name ## original name, e.g., ecoli
        self.url = url
        self.root = osp.join(root, name)
        self.kwargs = kwargs
        self.pdblist = os.path.join(os.path.dirname(__file__), '..', 'pkg_data', 'tm_pdblist.txt').readlines()
        self.n_prots = len(self.pdblist)

        super(PDBBindRefined, self).__init__(self.root, transform, pre_transform)

    @property
    def raw_file_names(self):
        if not osp.exists(self.raw_dir):
            return []
        return ['pdbs', 'tm_scores.pkl', 'rmsd.pkl']

    @property
    def processed_file_names(self):
        return [f'data_{i}.pt' for i in range(self.prots)]

    def download(self):
        print(f"Downloading to {self.root}...")
        tarball = wget.download(self.url, out=os.path.join(self.root))
        os.rename()

    def process(self):
        ### read pyg graph list
        """ Each graph will have an is_site tensor with 1 if the
        residue is in the binding site, 0 else.
        """
        data_list = []
        todo_pdbs = osp.join(self.raw_dir, 'pdbs')
        for i, pdb in tqdm(enumerate(todo_pdbs), total=len(todo_pdbs)):
            pdb_dir = osp.join(self.raw_dir, 'pdbs')
            graph = pdb2pyg(osp.join(pdb_dir, pdb), **self.kwargs)
            graph.name = pdb

            if self.pre_transform is not None:
                graph = self.pre_transform(graph)

            torch.save(graph, self.processed_paths[i])

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, f'data_{idx}.pt'))
        return data
