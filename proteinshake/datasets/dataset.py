# -*- coding: utf-8 -*-
"""
Base dataset class for protein 3D structures.
"""
import os, gzip, inspect, time, itertools, tarfile, io
from collections import defaultdict
from functools import cached_property
import multiprocessing as mp

import pandas as pd
import numpy as np
import freesasa
from biopandas.pdb import PandasPdb
from tqdm import tqdm
from joblib import Parallel, delayed
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph
from fastavro import reader as avro_reader

from proteinshake.transforms import IdentityTransform
from proteinshake.utils import download_url, save, load, unzip_file, write_avro, Generator

AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}
AA_ONE_TO_THREE = {v:k for k, v in AA_THREE_TO_ONE.items()}

# maps the date-format release to Zenodo identifier
RELEASES = {
    'latest': '1175739',
}

class Dataset():
    """ Base dataset class. Holds the logic for downloading and parsing PDB files.
    If ``use_precomputed=True``, fetched pre-processed data from Zenodo.
    Else, builds the dataset from scratch by executing: :meth:`download()` to fetch structures in PDB format, then :meth:`parse()` is applied to each to extract the relevant info and store it in a protein dictionary which has three outer keys ``'protein'``, ``'residue'``, and ``'atom'``. Subclassing :meth:`add_protein_attributes` lets the user include custom attributes.

    .. note::

        All child classes inherit these attributes and optionally add their own.

    .. list-table:: Annotations
      :widths: 25 35 45
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Protein identifier
        - :code:`protein['protein']['ID']`
        - ``'1JC8'``
      * - Sequence
        - :code:`protein['protein']['sequence']`
        - ``'IWGDSGKLITTA'``
      * - Assigned train/val/test split
        - :code:`protein['protein']['sequence_split_<CUTOFF>']`, :code:`protein['protein']['structure_split_<CUTOFF>']`
        - ``'train'``
      * - Residue position on chain
        - :code:`protein['residue']['residue_number']`
        - ``[2, 3, 4, 5, 6, 7, 8, 9, 10, ...,]``
      * - Amino acid type (single letter)
        - :code:`protein['residue']['residue_type']`
        - :code:`['I', 'W', 'G', 'D', 'S',..]`
      * - 3D coordinates
        - :code:`protein[{'residue' | 'atom'}][{'x'|'y'|'z'}]`
        - ``[5.191999912261963, 3.9860000610351562,..]``
      * - Solvent accessible surface area
        - :code:`protein[{'residue'|'atom'}]['SASA']`
        - :code:`[242.03138732910156, 136.46714782714844,... ]`
      * - Relative accessible surface area
        - :code:`protein['residue']['RSA']`
        - :code:`[1.377291202545166, 0.5476430058479309,... ]`
      * - Atom position
        - :code:`protein['atom']['atom_number']`
        - :code:`[7, 8, 9, 10, 11,..]`
      * - Atom type
        - :code:`protein['atom']['atom_type']`
        - :code:`['N', 'CA', 'C', 'O',...]`


    Arguments
    -----------
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
    minimum_length: int, default 10
        Proteins smaller than minimum_length residues will be skipped.
    maximum_length: int, default 2048
        Proteins larger than maximum_length residues will be skipped.
    exclude_ids: list, default []
        Exclude PDB IDs from the dataset.
    """

    additional_files = [] # indicates the additional file names that are to be included in the release
    exlude_args_from_signature = []

    def __init__(self,
            root                           = 'data',
            use_precomputed                = True,
            release                        = 'latest',
            only_single_chain              = False,
            check_sequence                 = False,
            n_jobs                         = 1,
            minimum_length                 = 10,
            maximum_length                 = 2048,
            exclude_ids                    = [],
            ):
        self.repository_url = f'https://sandbox.zenodo.org/record/{RELEASES[release]}/files'
        self.n_jobs = n_jobs
        self.use_precomputed = use_precomputed
        self.root = root
        self.minimum_length = minimum_length
        self.maximum_length = maximum_length
        self.only_single_chain = only_single_chain
        self.check_sequence = check_sequence
        self.release = release
        self.exclude_ids = exclude_ids
        
        os.makedirs(f'{self.root}', exist_ok=True)
        self.check_signature()

        if not use_precomputed:
            self.start_download()
            self.parse()
        else:
            self.check_signature_same_as_hosted()

    def compute_signature(self, use_defaults=False):
        signature = dict(inspect.signature(self.__init__).parameters.items())
        class_object = self.__class__
        while True: # add base signatures to subclass signature
            signature = {**dict(inspect.signature(class_object.__init__).parameters.items()), **signature}
            if len(class_object.__bases__) == 0: break
            class_object = class_object.__bases__[0]
        arg_names = [n for n in signature.keys() if not n in ['self', 'args', 'kwargs', 'n_jobs', 'root']+self.exlude_args_from_signature]
        if use_defaults:
            return self.name + ' | ' + ', '.join([k + '=' + str(signature[k].default) for k in arg_names])
        return self.name + ' | ' + ', '.join([k + '=' + str(getattr(self, k)) for k in arg_names])

    @cached_property
    def default_signature(self):
        return self.compute_signature(use_defaults=True)

    @cached_property
    def signature(self):
        return self.compute_signature(use_defaults=False)

    def check_signature(self):
        if os.path.exists(f'{self.root}/signature.txt'):
            with open(f'{self.root}/signature.txt','r') as file:
                assert file.read() == self.signature, 'The Dataset is called with different arguments than were used to create it. Delete or change the root.'
        else:
            with open(f'{self.root}/signature.txt','w') as file:
                file.write(self.signature)

    def check_signature_same_as_hosted(self):
        """ Safety check to ensure the provided dataset arguments are the same as were used to precompute the datasets. Only relevant with `use_precomputed=True`.
        """
        assert self.signature == self.default_signature, 'The dataset arguments do not match the precomputed dataset arguments (the default settings). Set use_precomputed to False if you wish to generate a new dataset.'

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


        .. code-block:: python

            >>> from proteinshake.datasets import RCSBDataset
            >>> protein = next(RCSBDataset().proteins())
        """
        self.download_precomputed(resolution=resolution)
        with open(f'{self.root}/{self.name}.{resolution}.avro', 'rb') as file:
            total = int(avro_reader(file).metadata['number_of_proteins'])
        def reader():
            with open(f'{self.root}/{self.name}.{resolution}.avro', 'rb') as file:
                for x in avro_reader(file):
                    yield x
        return Generator(reader(), total)

    @property
    def limit(self):
        """ Used only in testing, where this method is mock.patched to a small number. Default None.

        Returns
        -------
        int
            The limit to be applied to the number of downloaded/parsed files.
        """
        return None

    @property
    def name(self):
        return self.__class__.__name__

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
            A protein object. See :meth:`proteinshake.datasets.Dataset.parse_pdb` for details.

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
        if not os.path.exists(f'{self.root}/{self.name}.{resolution}.avro'):
            download_url(f'{self.repository_url}/{self.name}.{resolution}.avro.gz', f'{self.root}')
            print('Unzipping...')
            unzip_file(f'{self.root}/{self.name}.{resolution}.avro.gz')

    def parse(self):
        """ Parses all PDB files returned from :meth:`proteinshake.datasets.Dataset.get_raw_files()` and saves them to disk. Can run in parallel.
        """
        if os.path.exists(f'{self.root}/{self.name}.residue.avro'):
            return
        # parse and filter
        paths = self.get_raw_files()[:self.limit]
        proteins = Parallel(n_jobs=self.n_jobs)(delayed(self.parse_pdb)(path) for path in tqdm(paths, desc='Parsing'))
        before = len(proteins)
        proteins = [p for p in proteins if p is not None]
        print(f'Filtered {before-len(proteins)} proteins.')
        residue_proteins = [{'protein':p['protein'], 'residue':p['residue']} for p in proteins]
        atom_proteins = [{'protein':p['protein'], 'atom':p['atom']} for p in proteins]
        write_avro(residue_proteins, f'{self.root}/{self.name}.residue.avro')
        write_avro(atom_proteins, f'{self.root}/{self.name}.atom.avro')

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

        # add surface accessible area
        structure = freesasa.Structure(path)
        result = freesasa.calc(structure)
        residue_result = result.residueAreas()
        atom_sasa, residue_sasa, residue_rsa = [], [], []
        for i in atom_df['atom_number']:
            try:
                atom_sasa.append(result.atomArea(i))
            except:
                atom_sasa.append(-1)
        for i,chain in zip(residue_df['residue_number'], residue_df['chain_id']):
            try:
                residue_sasa.append(residue_result[chain][str(i)].total)
                residue_rsa.append(residue_result[chain][str(i)].relativeTotal)
            except:
                residue_sasa.append(-1)
                residue_rsa.append(-1)
            

        # create protein_dict
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
                'SASA': residue_sasa,
                'RSA': residue_rsa,
            },
            'atom': {
                'atom_number': atom_df['atom_number'].tolist(),
                'atom_type': atom_df['atom_type'].tolist(),
                'residue_number': atom_df['residue_number'].tolist(),
                'residue_type': atom_df['residue_type'].tolist(),
                'x': atom_df['x'].tolist(),
                'y': atom_df['y'].tolist(),
                'z': atom_df['z'].tolist(),
                'SASA': atom_sasa,
            },
        }

        # only include chains if multi-chain protein
        if not self.only_single_chain: 
            protein['residue']['chain_id'] = residue_df['chain_id'].tolist()
            protein['atom']['chain_id'] = atom_df['chain_id'].tolist()

        # add pLDDT from AlphaFold
        if self.name == 'AlphaFoldDataset':
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
        IONS = ['ZN', 'MG']
        df = PandasPdb().read_pdb_from_list(filtered_lines).df['ATOM']
        df = df.loc[~df['residue_name'].isin(IONS)]
        df['residue_name'] = df['residue_name'].map(lambda x: AA_THREE_TO_ONE[x] if x in AA_THREE_TO_ONE else None)
        #df['atom_name'] = df['atom_name'].map(lambda x: x[0]) # each atom is a multi-letter code where the first letter indicates the atom type
        df = df.rename(columns={
            'atom_name': 'atom_type',
            'residue_name': 'residue_type',
            'x_coord': 'x',
            'y_coord': 'y',
            'z_coord': 'z',
        })
        df = df.sort_values(by=['chain_id', 'residue_number', 'atom_number'])
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
        if len(df['residue_number'].unique()) > self.maximum_length:
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
    
    def to_graph(self, resolution='residue', transform=IdentityTransform(), *args, **kwargs):
        """ Converts the raw dataset to a graph dataset. See :meth:`proteinshake.representations.GraphDataset` for arguments.

        Returns
        -------
        proteinshake.representations.GraphDataset
            The dataset in graph representation.
        """
        from proteinshake.representations import GraphDataset
        proteins = self.proteins(resolution=resolution)
        return GraphDataset(Generator((transform(p) for p in proteins), len(proteins)),
                            self.root,
                            self.name,
                            resolution,
                            *args,
                            **kwargs)

    def to_point(self, resolution='residue', transform=IdentityTransform(), *args, **kwargs):
        """ Converts the raw dataset to a point cloud dataset. See :meth:`proteinshake.representations.PointDataset` for arguments.

        Returns
        -------
        proteinshake.representations.PointDataset
            The dataset in point cloud representation.
        """
        from proteinshake.representations import PointDataset
        proteins = self.proteins(resolution=resolution)
        return PointDataset(Generator((transform(p) for p in proteins), len(proteins)),
                            self.root,
                            self.name,
                            resolution,
                            *args,
                            **kwargs)

    def to_voxel(self, resolution='residue', transform=IdentityTransform(), *args, **kwargs):
        """ Converts the raw dataset to a voxel dataset. See :meth:`proteinshake.representations.VoxelDataset` for arguments.

        Returns
        -------
        proteinshake.representations.VoxelDataset
            The dataset in voxel representation.
        """
        from proteinshake.representations import VoxelDataset
        proteins = self.proteins(resolution=resolution)
        return VoxelDataset(Generator((transform(p) for p in proteins), len(proteins)),
                            self.root,
                            self.name,
                            resolution,
                            *args,
                            **kwargs)
