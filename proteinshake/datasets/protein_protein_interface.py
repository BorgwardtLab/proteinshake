# -*- coding: utf-8 -*-
import os
import glob
import pandas as pd
from biopandas.pdb import PandasPdb
from sklearn.neighbors import KDTree
import numpy as np

from proteinshake.datasets import Dataset
from proteinshake.utils import extract_tar, download_url

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


    def __init__(self, cutoff=6, version='2020', **kwargs):
        self.version = version
        self.cutoff = cutoff
        super().__init__(**kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')[:self.limit]

    def get_id_from_filename(self, filename):
        return filename[:4]

    def get_interfaces(self, protein, cutoff=6, resolution='residue'):
        """Obtain interfacing residues within a single structure of polymers. Uses
        KDTree data structure for vector search.

        Parameters
        ----------
        protein: dict
            Parsed protein dictionary.

        Returns
        --------
            `list`: indicator list for each residue with 0 if not in interface and 1 else.
        """

        #3-D KD tree
        df = pd.DataFrame({
            'x': protein[resolution]['x'],
            'y': protein[resolution]['y'],
            'z': protein[resolution]['z'],
            'chain': protein[resolution]['chain_id'],
            'residue_index':protein[resolution]['residue_number']
        })
        resi_df = df.groupby(['residue_index', 'chain']).mean().reset_index()
        resi_coords = np.array([resi_df['x'].tolist(), resi_df['y'].tolist(), resi_df['z'].tolist()]).T
        kdt = KDTree(resi_coords, leaf_size=1)

        query = kdt.query_radius(resi_coords, cutoff)
        interface = set()
        for i,result in enumerate(query):
            res_index = resi_df.iloc[i].name
            this_chain = resi_df.iloc[i].chain
            for r in result:
                that_resi = resi_df.iloc[r].name
                that_chain = resi_df.iloc[r].chain
                if this_chain != that_chain:
                    interface.add(res_index)
                    interface.add(that_resi)

        resi_interface = []
        for r in protein[resolution]['residue_number']:
            resi_interface.append(r in interface)
        return resi_interface

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw/files', extract_members=True)

    def add_protein_attributes(self, protein):
        interface_atoms = self.get_interfaces(protein, self.cutoff, resolution='atom')
        protein['atom']['is_interface'] = interface_atoms
        protein['residue']['is_interface'] = [
            val
            for val, atom_type in \
            zip(interface_atoms, protein['atom']['atom_type']) \
            if atom_type == 'CA'
        ]
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Protein-protein interface (residue-level)"
        desc['values'] = 2
        desc['type'] = "Binary"
        return desc
