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
from torch_geometric.data import Dataset
from biopandas.pdb import PandasPdb

from torch_pdb.utils.convert import pdb2pyg

class PDBBindRefined(Dataset):
    def __init__(self, name, root='/tmp/var', transform=None, pre_transform=None,
                 version='2020',
                 server_address="https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/",
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
        self.url = server_address + f"PDBbind_v{version}_refined.tar.gz"
        self.root = osp.join(root, name)
        self.kwargs = kwargs

        super(PDBBindRefined, self).__init__(self.root, transform, pre_transform)

        self.data, self.slices = torch.load(self.processed_paths[0])

    @property
    def raw_file_names(self):
        if not osp.exists(self.raw_dir):
            return []
        data_list_names = [
            pdb_file for pdb_file in os.listdir(self.raw_dir) if pdb_file not in ['index', 'readme']
            ]
        return data_list_names

    @property
    def processed_file_names(self):
        return [f'data_{i}.pt' for i in range(len(self.raw_file_names))]

    def download(self):
        print(f"Downloading to {self.root}...")
        tarball = wget.download(self.url, out=self.root)
        # tarball = download_url(self.url, self.root)
        print("Extracting...")
        extract_tar(osp.join(self.root, tarball), osp.join(self.root))
        print("Renaming")
        os.rename(osp.join(self.root, 'refined-set'), self.raw_dir)

    def process(self):
        ### read pyg graph list
        """ Each graph will have an is_site tensor with 1 if the
        residue is in the binding site, 0 else.
        """
        todo_pdbs = self.raw_file_names 
        data_list = []
        for i, pdb in tqdm(enumerate(todo_pdbs), total=len(todo_pdbs)):
            pdb_dir = osp.join(self.raw_dir, pdb)
            graph = pdb2pyg(osp.join(pdb_dir, f'{pdb}_protein.pdb'), **self.kwargs)
            graph.name = pdb

            is_site = torch.zeros_like(graph.residue_number)

            # get binding site residues
            struc = PandasPdb().read_pdb(osp.join(pdb_dir, f'{pdb}_pocket.pdb'))
            binding_residues = torch.tensor(struc.df['ATOM']['residue_number'].unique())
            is_site[(graph.residue_number.unsqueeze(1) == binding_residues).sum(dim=1).nonzero()] = 1.

            if self.pre_transform is not None:
                graph = self.pre_transform(graph)

            torch.save(graph, self.processed_paths[i])

    def len(self):
        return len(self.processed_file_names)

    def get(self, idx):
        data = torch.load(osp.join(self.processed_dir, f'data_{idx}.pt'))
        return data
