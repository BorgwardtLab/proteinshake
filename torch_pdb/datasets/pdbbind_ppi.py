# -*- coding: utf-8 -*-
import os
import glob

from torch_pdb.datasets import TorchPDBDataset
from torch_pdb.utils import get_interfaces, extract_tar, download_url

class PDBBindPPI(TorchPDBDataset):
    """Protein-protein complexes with annotated interfaces. Residues
    in each protein are marked with a boolean `is_interface` to indicate
    residues defined to belong to the interface of two protein chains.
    The default threshold for determining interface residues is 6 Angstroms.

    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.
    cutoff: float
        Distance in angstroms within which a pair of residues is considered to
        belong to the interface.
    """


    def __init__(self, cutoff=6, version='2020', **kwargs):
        self.version = version
        self.cutoff = cutoff
        super().__init__(**kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')[:self.download_limit()]

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        os.rename(f'{self.root}/raw/PP', f'{self.root}/raw/files')

    def add_protein_attributes(self, protein):
        protein['is_interface'] = get_interfaces(protein, self.cutoff)
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Protein-protein interface (residue-level)"
        desc['values'] = 2
        desc['type'] = "Binary"
        return desc
