# -*- coding: utf-8 -*-
"""
Base dataset class for protein 3D structures.
"""
import os, gzip, inspect

import pandas as pd
from biopandas.pdb import PandasPdb
from tqdm import tqdm
import numpy as np
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from joblib import Parallel, delayed

from torch_pdb.utils import download_url, save, load
from torch_pdb.representations import GraphDataset, PointDataset, VoxelDataset

three2one = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}


class TorchPDBDataset():
    """ Base dataset class. Holds the logic for downloading and parsing .pdb files.

    Parameters
    ----------
    root: str, default 'data'
        The data root directory to store both raw and parsed data.
    use_precomputed: bool, default True
        If `True`, will download the processed dataset from torch-pdb (recommended). If `False`, will download the raw data from the original sources and process them on your device. You can use this option if you wish to create a custom dataset. Using `False` is compute-intensive, consider increasing `n_jobs`.
    release: str, default '12JUL2022'
        The GitHub tag of the dataset release. See https://github.com/BorgwardtLab/torch-pdb/releases for all available releases. Latest (default) is recommended.
    only_single_chain: bool, default False
        If `True`, will only use single-chain proteins.
    check_sequence: bool, default False
        If `True`, will discard proteins whose primary sequence is not identical with the sequence of amino acids in the structure. This can happen if the structure is not complete (e.g. for parts that could not be crystallized).
    n_jobs: int, default 1
        The number of jobs for downloading and parsing files. It is recommended to increase the number of jobs with `use_precomputed=False`.
    """
    def __init__(self,
            root                = 'data',
            use_precomputed     = True,
            release             = '12JUL2022'
            only_single_chain   = False,
            check_sequence      = False,
            n_jobs              = 1,
            ):
        self.n_jobs = n_jobs
        self.use_precomputed = use_precomputed
        self.root = root
        self.only_single_chain = only_single_chain
        self.check_sequence = check_sequence
        self.release = release
        self.check_arguments_same_as_hosted()
        self._download()
        self.parse()
        self.proteins = load(f'{self.root}/{self.__class__.__name__}.json.gz')

    def download_limit(self):
        """ Used only in testing, where this method is mock.patched to a small number.

        Returns
        -------
        None
            The limit to be applied to the number of downloaded/parsed files.
        """
        return None

    def check_arguments_same_as_hosted(self):
        """ Safety check to ensure the provided dataset arguments are the same as were used to precompute the datasets. Only relevant with `use_precomputed=True`.
        """
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
        """ Implement me in a subclass!

        Returns a list of all valid PDB file paths for this dataset. Usually takes the form `glob.glob(f'{self.root}/raw/files/*.pdb')` to search for all pdb files in the root, but can be different in some cases (e.g. with pdb.gz files).

        Returns
        -------
        list
            The list of raw PDB files used in this dataset.
        """
        raise NotImplementedError

    def get_id_from_filename(self, filename):
        """ Implement me in a subclass!

        Extracts an identifier from the pdb filename. This identifier is used in the `ID` field of the parsed protein object. Usually something like `filename[:4]`.

        Parameters
        ----------
        filename: str
            Path to a PDB file.

        Returns
        -------
        str
            A PDB identifier or other ID.
        """
        raise NotImplementedError

    def download(self):
        """ Implement me in a subclass!

        This method is responsible for downloading and extracting raw pdb files from a databank source. All PDB files should be dumped in `f'{self.root}/raw/files`. See e.g. PDBBindRefined for an example.
        """
        raise NotImplementedError

    def add_protein_attributes(self, protein):
        """ Implement me in a subclass!

        This method annotates protein objects with addititional information, such as functional labels or classes. It takes a protein object (a dictionary), modifies, and returns it. Usually, this would utilize the `ID` attribute to load an annotation file or to query information from a database.

        Parameters
        ----------
        protein: dict
            A protein object with `ID`, `sequence`, `coords` and other features. See `TorchPDBDataset.parse_pdb()` for details.

        Returns
        -------
        dict
            The protein object with a new attribute added.
        """
        return protein

    def download_complete(self):
        """ Dumps a marker file when the download was successful, to skip downloading next time.
        """
        with open(f'{self.root}/raw/done.txt','w') as file:
            file.write('done.')

    def _download(self):
        """ Helper function to prepare the download. Switches between precomputed and raw download. Creates necessary subdirectories.
        """
        if self.use_precomputed:
            self.download_precomputed()
        else:
            if os.path.exists(f'{self.root}/raw/done.txt'):
                return
            os.makedirs(f'{self.root}/raw/files', exist_ok=True)
            self.download()
            self.download_complete()

    def download_precomputed(self):
        """ Downloads the precomputed dataset from torch-pdb.
        """
        os.makedirs(f'{self.root}', exist_ok=True)
        if not os.path.exists(f'{self.root}/{self.__class__.__name__}.json.gz'):
            download_url(f'https://github.com/BorgwardtLab/torch-pdb/releases/download/{self.release}/{self.__class__.__name__}.json.gz', f'{self.root}')

    def parse(self):
        """ Parses all PDB files returned from `self.get_raw_files()` and saves them to disk. Can run in parallel.
        """
        if os.path.exists(f'{self.root}/{self.__class__.__name__}.json.gz'):
            return
        proteins = Parallel(n_jobs=self.n_jobs)(delayed(self.parse_pdb)(path) for path in tqdm(self.get_raw_files(), desc='Parsing PDB files'))
        #proteins = [self.parse_pdb(path) for path in tqdm(self.get_raw_files(), desc='Parsing PDB files')]
        before = len(proteins)
        proteins = [p for p in proteins if p is not None]
        print(f'Filtered {before-len(proteins)} proteins.')
        save(proteins, f'{self.root}/{self.__class__.__name__}.json.gz')

    def parse_pdb(self, path):
        """ Parses a single PDB file first into a DataFrame, then into a protein object (a dictionary). Also validates the PDB file and provides the hook for `add_protein_attributes`. Should return `None` if the protein was found to be invalid.

        Parameters
        ----------
        path: str
            Path to PDB file.

        Returns
        -------
        dict
            A protein object.
        """
        df = self.pdb2df(path)
        if not self.validate(df):
            return None
        protein = {
            'ID': self.get_id_from_filename(os.path.basename(path)),
            'sequence': ''.join(df['residue_name']),
            'residue_index': df['residue_number'].tolist(),
            'coords': df[['x_coord','y_coord','z_coord']].values.tolist(),
        }
        if not self.only_single_chain: # only include chains if multi-chain protein
            protein['chain_id'] = df['chain_id'].tolist()
        protein = self.add_protein_attributes(protein)
        return protein

    def pdb2df(self, path):
        """ Parses a single PDB file to a DataFrame (with biopandas). Also deals with multiple structure models in a PDB (e.g. from NMR) by only selecting the first model.

        Parameters
        ----------
        path: str
            Path to PDB file.

        Returns
        -------
        DataFrame
            A biopandas DataFrame of the PDB file.
        """
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
        """ Performs several sanity checks for validity of the protein DataFrame.

        Parameters
        ----------
        df: DataFrame
            A protein DataFrame.

        Returns
        -------
        bool
            Whether or not the DataFrame is valid.
        """
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

    def describe(self):
        """ Produces dataset statistics.

        Returns
        -------
        dict
            A dictionary of summary statistics of this dataset.
        """
        n_resi = len(self.data.residue_index) / len(self.data.ID)
        data = {'name': type(self).__name__,
                'num_proteins': len(self),
                'avg size (# residues)': n_resi
               }
        return data

    def to_graph(self, *args, **kwargs):
        """ Converts the raw dataset to a graph dataset. See `GraphDataset` for arguments.

        Returns
        -------
        GraphDataset
            The dataset in graph representation.
        """
        return GraphDataset(self.root, self.proteins, *args, **kwargs)

    def to_point(self, *args, **kwargs):
        """ Converts the raw dataset to a point cloud dataset. See `PointDataset` for arguments.

        Returns
        -------
        PointDataset
            The dataset in point cloud representation.
        """
        return PointDataset(self.root, self.proteins, *args, **kwargs)

    def to_voxel(self, *args, **kwargs):
        """ Converts the raw dataset to a voxel dataset. See `VoxelDataset` for arguments.

        Returns
        -------
        VoxelDataset
            The dataset in voxel representation.
        """
        return VoxelDataset(self.root, self.proteins, *args, **kwargs)
