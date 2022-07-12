"""
Base dataset class for 3D structures.
"""
# -*- coding: utf-8 -*-
import os, gzip, inspect

import torch
import pandas as pd
from biopandas.pdb import PandasPdb
from tqdm import tqdm
import numpy as np
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from joblib import Parallel, delayed
from torch_geometric.data import InMemoryDataset, Data, extract_tar, download_url
from torch_geometric.utils import from_scipy_sparse_matrix

from torch_pdb.utils import one_hot

three2one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

class TorchPDBDataset(InMemoryDataset):
    """ Base dataset class for holding 3D structures and encoding as graphs."""
    def __init__(self,
            root                = 'data',
            name                = 'proteins',
            node_embedding      = one_hot,
            graph_construction  = 'eps',
            eps                 = 8,
            k                   = 5,
            weighted_edges      = False,
            only_single_chain   = False,
            check_sequence      = False,
            n_jobs              = 1,
            use_precomputed     = True,
            release             = 'JUL2022'
            ):
        self.n_jobs = n_jobs
        self.use_precomputed = use_precomputed
        self.root = root
        self.name = name
        self.node_embedding = node_embedding
        self.graph_construction = graph_construction
        self.eps = eps
        self.k = k
        self.weighted_edges = weighted_edges
        self.only_single_chain = only_single_chain
        self.check_sequence = check_sequence
        self.release = release
        self.check_arguments_same_as_hosted()
        super().__init__(root)
        self._download() # some weird quirk requires this if .download() / .process() is not defined on the lowest inheritance level, might want to look into this at some point
        self._process()
        self.data, self.slices = torch.load(f'{self.root}/processed/{self.name}.pt')

    def check_arguments_same_as_hosted(self):
        signature = inspect.signature(self.__init__)
        default_args = {
            k: v.default
            for k, v in signature.parameters.items()
            if v.default is not inspect.Parameter.empty
            and (self.__class__.__name__ != 'AlphaFoldDataset' or k != 'organism')
        }
        if self.__class__.__bases__[0].__name__ != 'TorchPDBDataset':
            signature = inspect.signature(self.__class__.__bases__[0].__init__)
            super_args = {
                k: v.default
                for k, v in signature.parameters.items()
                if v.default is not inspect.Parameter.empty
            }
            default_args = {**super_args, **default_args}
        if self.use_precomputed and not all([v == getattr(self, k) for k,v in default_args.items()]):
            print('Error: The dataset arguments do not match the precomputed dataset arguments (the default settings). Set use_precomputed to False if you wish to generate a new dataset.')
            exit()

    def get_raw_files(self):
        ''' Returns a list of all valid PDB files.
        Implement me! '''
        raise NotImplementedError

    def get_id_from_filename(self, filename):
        ''' Takes in raw filename `xyz_abc.pdb` and returns a PDBID.
        Implement me! '''
        raise NotImplementedError

    def download(self):
        ''' Dumps data to /raw and /raw/files/*.pdb.
        Implement me! '''
        raise NotImplementedError

    def describe(self):
        """ Produce dataset statistics.
        """
        n_resi = len(self.data.residue_index) / len(self.data.ID)
        data = {'name': type(self).__name__,
                'num_proteins': len(self),
                'avg size (# residues)': n_resi
               }
        return data

    def add_protein_attributes(self, protein):
        ''' Implement me! '''
        return protein

    def download_complete(self):
        print('Download complete.')
        with open(f'{self.root}/raw/done.txt','w') as file:
            file.write('done.')

    def _process(self):
        if os.path.exists(f'{self.root}/processed/{self.name}.pt'):
            return
        if not os.path.exists(f'{self.root}/{self.__class__.__name__}.pt'):
            self.parse()
        os.makedirs(f'{self.root}/processed', exist_ok=True)
        self.process()

    def _download(self):
        if self.use_precomputed:
            self.download_precomputed()
        else:
            if os.path.exists(f'{self.root}/raw/done.txt'):
                return
            os.makedirs(f'{self.root}/raw/files', exist_ok=True)
            self.download()
            self.download_complete()
            self.parse()

    def download_precomputed(self):
        if not os.path.exists(f'{self.root}/{self.__class__.__name__}.pt'):
            download_url(f'https://github.com/BorgwardtLab/torch-pdb/releases/download/{self.release}/{self.__class__.__name__}.pt', f'{self.root}')

    def parse(self):
        #proteins = Parallel(n_jobs=self.n_jobs)(delayed(self.parse_pdb)(path) for path in tqdm(self.get_raw_files(), desc='Parsing PDB files'))
        proteins = [self.parse_pdb(path) for path in tqdm(self.get_raw_files(), desc='Parsing PDB files')]
        proteins = [p for p in proteins if p is not None]
        torch.save(proteins, f'{self.root}/{self.__class__.__name__}.pt')

    def process(self):
        proteins = torch.load(f'{self.root}/{self.__class__.__name__}.pt')
        convert = lambda p: self.graph2pyg(self.protein2graph(p), info=p)
        data_list = Parallel(n_jobs=self.n_jobs)(delayed(convert)(p) for p in tqdm(proteins, desc='Converting proteins to graphs'))
        print('Saving...')
        data, slices = self.collate(data_list)
        torch.save((data, slices), f'{self.root}/processed/{self.name}.pt')
        print('Dataset ready.')

    def parse_pdb(self, path):
        df = self.pdb2df(path)
        if not self.validate(df):
            return None
        protein = {
            'ID': self.get_id_from_filename(os.path.basename(path)),
            'sequence': ''.join(df['residue_name']),
            'residue_index': torch.tensor(df['residue_number'].tolist()).int(),
            'chain_id': df['chain_id'].tolist(),
            'coords': torch.stack([
                torch.tensor(df['x_coord'].to_list()),
                torch.tensor(df['y_coord'].to_list()),
                torch.tensor(df['z_coord'].to_list())
            ], dim=1).long(),
        }
        protein = self.add_protein_attributes(protein)
        return protein

    def pdb2df(self, path):
        if path.endswith('.gz'):
            with gzip.open(path, 'rb') as file:
                lines = file.read().decode('utf-8').split('\n')
        else:
            with open(path, 'r') as file:
                lines = file.read().split('\n')
        # filter only the first model
        filtered_lines, in_model, model_done = [], False, False
        for line in lines:
            if line.startswith('MODEL'):
                in_model = True
            if in_model and model_done:
                continue
            if line.startswith('ENDMDL'):
                model_done = True
                in_model = False
            filtered_lines.append(line)
        df = PandasPdb().read_pdb_from_list(filtered_lines).df['ATOM']
        df = df[df['atom_name'] == 'CA']
        df['residue_name'] = df['residue_name'].map(lambda x: three2one[x] if x in three2one else None)
        df = df.sort_values('residue_number')
        return df

    def validate(self, df):
        # check if single chain protein
        if self.only_single_chain and len(df['chain_id'].unique()) > 1:
            return False
        # check if sequence and structure are consistent
        if self.check_sequence and not np.array_equal(df.index, np.arange(1,len(df)+1)):
            return False
        # check if all standard amino acids
        if not all(df['residue_name'].map(lambda x: not x is None)):
            return False
        # check if sequence is empty (e.g. with non-canonical amino acids)
        if not sum(df['residue_name'].map(lambda x: not x is None)) > 0:
            return False
        return True

    def protein2graph(self, protein):
        nodes = self.node_embedding(protein['sequence'])
        if self.graph_construction == 'eps':
            mode = 'distance' if self.weighted_edges else 'connectivity'
            adj = radius_neighbors_graph(protein['coords'], radius=self.eps, mode=mode)
        elif self.graph_construction == 'knn':
            adj = kneighbors_graph(protein['coords'], k=self.k)
        return (nodes, adj)

    def graph2pyg(self, graph, info={}):
        nodes = torch.Tensor(graph[0]).float()
        edges = from_scipy_sparse_matrix(graph[1])
        return Data(x=nodes, edge_index=edges[0].long(), edge_attr=edges[1].unsqueeze(1).float(), **info)
