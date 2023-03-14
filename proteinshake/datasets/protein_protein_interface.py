# -*- coding: utf-8 -*-
import os
import glob

from proteinshake.datasets import Dataset
from proteinshake.utils import get_interfaces, extract_tar, download_url

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

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        os.rename(f'{self.root}/raw/PP', f'{self.root}/raw/files')

    def add_protein_attributes(self, protein):
        interface_atoms = get_interfaces(protein, self.cutoff, resolution='atom')
        protein['atom']['is_interface'] = interface_atoms
        protein['residue']['is_interface'] = [val
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
