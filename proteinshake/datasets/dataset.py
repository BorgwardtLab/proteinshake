# -*- coding: utf-8 -*-
"""
Base dataset class for protein 3D structures.
"""
import os, gzip, inspect, time, itertools, tarfile, io

import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
from tqdm import tqdm
from joblib import delayed
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from fastavro import reader as avro_reader

from proteinshake.utils import (download_url,
                                save,
                                load,
                                unzip_file,
                                ProgressParallel,
                                write_avro,
                                )

AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
AA_ONE_TO_THREE = {v:k for k, v in AA_THREE_TO_ONE.items()}

# maps the date-format release to Zenodo identifier
RELEASES = {
    'latest': '1132621',
    '23SEP2022': '1132621'
}

class Dataset():
    """ Base dataset class. Holds the logic for downloading and parsing PDB files.

    Parameters
    ----------
    root: str, default 'data'
        The data root directory to store both raw and parsed data.
    use_precomputed: bool, default True
        If `True`, will download the processed dataset from the ProteinShake repository (recommended). If `False`, will download the raw data from the original sources and process them on your device. You can use this option if you wish to create a custom dataset. Using `False` is compute-intensive, consider increasing `n_jobs`.
    release: str, default '12JUL2022'
        The tag of the dataset release. See https://github.com/BorgwardtLab/proteinshake/releases for all available releases. "latest" (default) is recommended.
    only_single_chain: bool, default False
        If `True`, will only use single-chain proteins.
    check_sequence: bool, default False
        If `True`, will discard proteins whose primary sequence is not identical with the sequence of amino acids in the structure. This can happen if the structure is not complete (e.g. for parts that could not be crystallized).
    n_jobs: int, default 1
        The number of jobs for downloading and parsing files. It is recommended to increase the number of jobs with `use_precomputed=False`.
    min_size: int, default 10
        Proteins smaller than min_size residues will be skipped.
    """
    def __init__(self,
            root                = 'data',
            use_precomputed     = True,
            release             = 'latest',
            only_single_chain   = False,
            check_sequence      = False,
            n_jobs              = 1,
            minimum_length      = 10,
            exclude_ids         = []
            ):
        self.repository_url = f'https://sandbox.zenodo.org/record/{RELEASES[release]}/files'
        self.n_jobs = n_jobs
        self.use_precomputed = use_precomputed
        self.root = root
        self.minimum_length = minimum_length
        self.only_single_chain = only_single_chain
        self.check_sequence = check_sequence
        self.release = release
        self.exclude_ids = exclude_ids
        os.makedirs(f'{self.root}', exist_ok=True)
        if not use_precomputed:
            self.start_download()
            self.parse()
        else:
            self.check_arguments_same_as_hosted()

    def proteins(self, resolution='residue'):
        """ Returns a generator of proteins from the avro file.

        Parameters
        ----------
        resolution: str, default 'residue'
            The resolution of the proteins. Can be 'atom' or 'residue'.

        Returns
        -------
        generator
            An avro reader object.

        int
            The total number of proteins in the file.
        """
        self.download_precomputed(resolution=resolution)
        with open(f'{self.root}/{self.__class__.__name__}.{resolution}.avro', 'rb') as file:
            total = int(avro_reader(file).metadata['number_of_proteins'])
        def reader():
            with open(f'{self.root}/{self.__class__.__name__}.{resolution}.avro', 'rb') as file:
                for x in avro_reader(file):
                    yield x
        return reader(), total

    def download_limit(self):
        """ Used only in testing, where this method is mock.patched to a small number. Default None.

        Returns
        -------
        int
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
            and (self.__class__.__name__ != 'Atom3DDataset' or k != 'atom_dataset')
        }
        if self.__class__.__bases__[0].__name__ != 'Dataset':
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
            A protein object. See `Dataset.parse_pdb()` for details.

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

    def start_download(self):
        """ Helper function to prepare the download. Creates necessary subdirectories.
        """
        if os.path.exists(f'{self.root}/raw/done.txt'):
            return
        os.makedirs(f'{self.root}/raw/files', exist_ok=True)
        self.download()
        self.download_complete()

    def download_precomputed(self, resolution='residue'):
        """ Downloads the precomputed dataset from the ProteinShake repository.
        """
        if not os.path.exists(f'{self.root}/{self.__class__.__name__}.{resolution}.avro'):
            download_url(f'{self.repository_url}/{self.__class__.__name__}.{resolution}.avro.gz', f'{self.root}')
            print('Unzipping...')
            unzip_file(f'{self.root}/{self.__class__.__name__}.{resolution}.avro.gz')

    def parse(self):
        """ Parses all PDB files returned from `self.get_raw_files()` and saves them to disk. Can run in parallel.
        """
        if os.path.exists(f'{self.root}/{self.__class__.__name__}.residue.avro'):
            return

        # parse and filter
        paths = self.get_raw_files()
        chunk_size = 1000
        chunks = [paths[i:i+chunk_size] for i in range(0,len(paths), chunk_size)]
        def parse_pdbs(chunk):
            return [self.parse_pdb(path) for path in chunk]
        n_jobs = 1#min(self.n_jobs, len(chunks)) # Parallelization does not really work for some reason
        if n_jobs == 1:
            proteins = [self.parse_pdb(path) for path in tqdm(paths, desc='Parsing')]
        else:
            proteins = ProgressParallel(n_jobs=n_jobs, total=len(chunks), desc='Parsing')(delayed(parse_pdbs)(chunk) for chunk in chunks)
            proteins = list(itertools.chain(*proteins))
        before = len(proteins)
        proteins = [p for p in proteins if p is not None]
        print(f'Filtered {before-len(proteins)} proteins.')
        residue_proteins = [{'protein':p['protein'], 'residue':p['residue']} for p in proteins]
        atom_proteins = [{'protein':p['protein'], 'atom':p['atom']} for p in proteins]
        write_avro(residue_proteins, f'{self.root}/{self.__class__.__name__}.residue.avro')
        write_avro(atom_proteins, f'{self.root}/{self.__class__.__name__}.atom.avro')

    def parse_pdb(self, path):
        """ Parses a single PDB file first into a DataFrame, then into a protein object (a dictionary). Also validates the PDB file and provides the hook for `add_protein_attributes`. Returns `None` if the protein was found to be invalid.

        Parameters
        ----------
        path: str
            Path to PDB file.

        Returns
        -------
        dict
            A protein object.
        """
        pdbid = self.get_id_from_filename(os.path.basename(path))
        if pdbid in self.exclude_ids:
            return None
        atom_df = self.pdb2df(path)
        residue_df = atom_df[atom_df['atom_type'] == 'CA']
        if not self.validate(atom_df):
            return None
        protein = {
            'protein': {
                'ID': pdbid,
                'sequence': ''.join(residue_df['residue_type']),
            },
            'residue': {
                'residue_number': residue_df['residue_number'].tolist(),
                'residue_type': residue_df['residue_type'].tolist(),
                'x': residue_df['x'].tolist(),
                'y': residue_df['y'].tolist(),
                'z': residue_df['z'].tolist(),
            },
            'atom': {
                'atom_number': atom_df['atom_number'].tolist(),
                'atom_type': atom_df['atom_type'].tolist(),
                'residue_number': atom_df['residue_number'].tolist(),
                'residue_type': atom_df['residue_type'].tolist(),
                'x': atom_df['x'].tolist(),
                'y': atom_df['y'].tolist(),
                'z': atom_df['z'].tolist(),
            },
        }
        if not self.only_single_chain: # only include chains if multi-chain protein
            protein['residue']['chain_id'] = residue_df['chain_id'].tolist()
            protein['atom']['chain_id'] = atom_df['chain_id'].tolist()
        # add pLDDT from AlphaFold
        if self.__class__.__name__ == 'AlphaFoldDataset':
            protein['residue']['pLDDT'] = residue_df['b_factor'].tolist()
            protein['atom']['pLDDT'] = atom_df['b_factor'].tolist()
        # add attributes
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
        df['residue_name'] = df['residue_name'].map(lambda x: AA_THREE_TO_ONE[x] if x in AA_THREE_TO_ONE else None)
        #df['atom_name'] = df['atom_name'].map(lambda x: x[0]) # each atom is a multi-letter code where the first letter indicates the atom type
        df = df.sort_values('atom_number')
        df = df.rename(columns={
            'atom_name': 'atom_type',
            'residue_name': 'residue_type',
            'x_coord': 'x',
            'y_coord': 'y',
            'z_coord': 'z',
        })
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
        if len(df['residue_number'].unique()) < self.minimum_length:
            return False
        # check if single chain protein
        if self.only_single_chain and len(df['chain_id'].unique()) > 1:
            return False
        # check if sequence and structure are consistent
        if self.check_sequence and not max(df['residue_number']) == len(df['residue_number'].unique()):
            return False
        # check if all standard amino acids
        if not all(df['residue_type'].map(lambda x: not x is None)):
            return False
        # check if sequence is empty (e.g. with non-canonical amino acids)
        if not sum(df['residue_type'].map(lambda x: not x is None)) > 0:
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

    def to_graph(self, resolution='residue', *args, **kwargs):
        """ Converts the raw dataset to a graph dataset. See `GraphDataset` for arguments.

        Returns
        -------
        GraphDataset
            The dataset in graph representation.
        """
        from proteinshake.representations import GraphDataset
        return GraphDataset(*self.proteins(resolution), self.root, resolution, *args, **kwargs)

    def to_point(self, resolution='residue', *args, **kwargs):
        """ Converts the raw dataset to a point cloud dataset. See `PointDataset` for arguments.

        Returns
        -------
        PointDataset
            The dataset in point cloud representation.
        """
        from proteinshake.representations import PointDataset
        return PointDataset(*self.proteins(resolution), self.root, resolution, *args, **kwargs)

    def to_voxel(self, resolution='residue', *args, **kwargs):
        """ Converts the raw dataset to a voxel dataset. See `VoxelDataset` for arguments.

        Returns
        -------
        VoxelDataset
            The dataset in voxel representation.
        """
        from proteinshake.representations import VoxelDataset
        return VoxelDataset(*self.proteins(resolution), self.root, resolution, *args, **kwargs)

    def to_surface(self, resolution='residue', *args, **kwargs):
        """ Converts the raw dataset to a voxel dataset. See `VoxelDataset` for arguments.

        Returns
        -------
        SurfaceDataset
            The dataset in surface representation.
        """
        from proteinshake.representations import SurfaceDataset
        return SurfaceDataset(*self.proteins(resolution), self.root, resolution, *args, **kwargs)
