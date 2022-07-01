# -*- coding: utf-8 -*-
import torch, os
from torch_geometric.data import InMemoryDataset, Data
from torch_geometric.utils import from_scipy_sparse_matrix
from biopandas.pdb import PandasPdb
from tqdm import tqdm
import numpy as np
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from torch_pdb.embeddings import one_hot
from joblib import Parallel, delayed

three2one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

class TorchPDBDataset(InMemoryDataset):
    def __init__(self,
            root,
            name,
            node_embedding      = one_hot,
            graph_construction  = 'eps',
            eps                 = 8,
            k                   = 5,
            weighted_edges      = False,
            only_single_chain   = False,
            check_sequence      = False,
            n_jobs              = 1,
            use_precomputed     = True,
            ):
        if not use_precomputed and n_jobs == 1:
            print('Downloading and processing an entire dataset with use_precompute = False is very slow. Consider increasing n_jobs.')
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
        super().__init__(root)
        self._download() # some weird quirk requires this if .download() / .process() is not defined on the lowest inheritance level, might want to look into this at some point
        self._process()
        self.data, self.slices = torch.load(self.processed_paths[0])

    def get_raw_files(self):
        ''' Implement me! '''
        raise NotImplementedError

    def get_id_from_filename(self, filename):
        ''' Implement me! '''
        raise NotImplementedError

    def download(self):
        ''' Implement me! '''
        raise NotImplementedError

    def add_protein_attributes(self, protein):
        ''' Implement me! '''
        return protein

    @property
    def raw_file_names(self):
        return ['done.txt']

    @property
    def processed_file_names(self):
        return [f'{self.name}.pt']

    def download_complete(self):
        print('Download complete.')
        with open(f'{self.root}/raw/{self.raw_file_names[0]}','w') as file:
            file.write('done.')

    def process(self):
        proteins = Parallel(n_jobs=self.n_jobs)(delayed(self.parse_pdb)(path) for path in tqdm(self.get_raw_files(), desc='Parsing PDB files'))
        proteins = [p for p in proteins if p is not None]
        convert = lambda p: self.graph2pyg(self.protein2graph(p), info=p)
        #data_list = Parallel(n_jobs=self.n_jobs)(delayed(convert)(p) for p in tqdm(proteins, desc='Converting proteins to graphs'))
        data_list = [convert(p) for p in tqdm(proteins)]
        print('Saving...')
        data, slices = self.collate(data_list)
        torch.save((data, slices), self.processed_paths[0])
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
        df = PandasPdb().read_pdb(path).df['ATOM']
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
