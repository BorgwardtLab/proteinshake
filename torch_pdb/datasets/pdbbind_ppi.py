# -*- coding: utf-8 -*-
import os
import glob

import torch
from torch_geometric.data import extract_tar, download_url

from torch_pdb.datasets import TorchPDBDataset
from torch_pdb.utils import get_interfaces

class PDBBindPPI(TorchPDBDataset):
    """"Dataset conatining proteins bound to small molecules. Residues
    in each protein are marked with a boolean `binding_site` to indicate
    residues defined to be inside the binding pocket.

    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.
    """


    def __init__(self, version='2020', **kwargs):
        self.version = version
        super().__init__(**kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')
    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        os.rename(f'{self.root}/raw/PP', f'{self.root}/raw/files')

    def add_protein_attributes(self, protein):
        protein['is_interface'] = get_interfaces(protein)
        return protein
