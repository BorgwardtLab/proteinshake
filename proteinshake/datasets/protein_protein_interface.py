# -*- coding: utf-8 -*-
import os
import glob
from joblib import Parallel, delayed
from collections import defaultdict, Counter

import pandas as pd
from biopandas.pdb import PandasPdb
from sklearn.neighbors import KDTree
import numpy as np

from proteinshake.datasets import Dataset
from proteinshake.utils import extract_tar, download_url, progressbar, load, save

class ProteinProteinInterfaceDataset(Dataset):
    """Protein-protein complexes from PDBBind with annotated interfaces. Residues
    and atoms in each protein are marked with a boolean `is_interface` to indicate
    residues/atoms defined to belong to the interface of two protein chains.
    The default threshold for determining interface residues is 6 Angstroms.
    See :meth:`proteinshake.utils.get_interfaces` for details.

    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.
    version: str
        PDBBind database version to use.
    cutoff: float
        Distance in angstroms within which a pair of residues is considered to
        belong to the interface.

    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 2839


   .. list-table:: Annotations
      :widths: 25 55 20
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Chain Identifier 
        - :code:`protein[{'residue' | 'atom'}]['chain_id']`
        - ``['X', 'X', ..'Y', 'Y']``
      * - Binding interface
        - :code:`protein[{'residue' | 'atom'}]['is_interface']`
        - ``[0, 0, .., 1, 0]``


    """
    additional_files = [
        'ProteinProteinInterfaceDataset.interfaces.json',
    ]

    def __init__(self, cutoff=6, version='2020', split_chains=True, **kwargs):
        self.version = version
        self.cutoff = cutoff
        kwargs['split_chains'] = split_chains
        super().__init__(**kwargs)

        def download_file(filename):
            if not os.path.exists(f'{self.root}/{filename}'):
                download_url(f'{self.repository_url}/{filename}.gz', f'{self.root}', verbosity=0)
                unzip_file(f'{self.root}/{filename}.gz')
            return load(f'{self.root}/{filename}')


        if not self.use_precomputed:
            self.parse_interfaces()

        self._interfaces = download_file(f'{self.name}.interfaces.json')

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/PP/*.pdb')[:self.limit]

    def get_id_from_filename(self, filename):
        return filename[:4]

    def get_contacts(self, protein, cutoff=6):
        """Obtain interfacing residues within a single structure of polymers. Uses
        KDTree data structure for vector search.

        Parameters
        ----------
        protein: dict
            Parsed protein dictionary.

        Returns
        --------
            `dict`: 2-level dictionary mapping a pair of chains to the list of interfacing residue positions (e.g `interfaces['A']['B'] = {(1, 3), (2, 4)}`
                    says that residues 1 and 2 of chain A are in contact with 3 and 4 in chain B. The positions are _indices_ in the residue/atom list.
        """

        def get_coords(p):
            return np.array([p['residue']['x'],
                             p['residue']['y'],
                             p['residue']['z']]).T


        def defaultdict_to_dict(default_dict):
            if isinstance(default_dict, defaultdict):
                default_dict = {k: defaultdict_to_dict(v) for k, v in default_dict.items()}
            elif isinstance(default_dict, dict):
                default_dict = {k: defaultdict_to_dict(v) for k, v in default_dict.items()}
            return default_dict

        interfaces = defaultdict(lambda: defaultdict(list))
        coords = get_coords(protein)
        kdt = KDTree(coords, leaf_size=1)

        seq_inds = [0]
        current_chain = protein['residue']['chain_id'][0]

        # get a vector of sequence positions
        ind = 0
        for chain_id in protein['residue']['chain_id'][1:]:
            if chain_id != current_chain:
                ind = -1
                current_chain = chain_id
            ind += 1
            seq_inds.append(ind)

        # if protein['protein']['ID'] == '2xns':
        print(protein['protein']['ID'], Counter(protein['residue']['chain_id']))

        query = kdt.query_radius(coords, cutoff)
        interface = []
        for i,result in enumerate(query):
            this_chain = protein['residue']['chain_id'][i]
            this_pos = seq_inds[i]
            for r in result:
                that_chain = protein['residue']['chain_id'][r]
                that_pos = seq_inds[r]
                if this_chain != that_chain:
                    # ugly , I know
                    interfaces[this_chain][that_chain].append((this_pos, that_pos))
                    interfaces[that_chain][this_chain].append((that_pos, this_pos))

        return defaultdict_to_dict(interfaces)

    def parse_interfaces(self):
        """ Get all interfaces and store in dict"""
        protein_dfs = Parallel(n_jobs=self.n_jobs)(delayed(self.parse_pdb)(path) for path in progressbar(self.get_raw_files(), desc='Parsing'))
        interfaces = {p['protein']['ID']: self.get_contacts(p, cutoff=self.cutoff) for p in protein_dfs}
        save(interfaces, f'{self.root}/{self.name}.interfaces.json')

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw/files', extract_members=True)

    def describe(self):
        desc = super().describe()
        desc['property'] = "Protein-protein interface (residue-level)"
        desc['values'] = 2
        desc['type'] = "Binary"
        return desc
