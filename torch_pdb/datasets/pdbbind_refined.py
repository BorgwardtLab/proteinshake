# -*- coding: utf-8 -*-
import glob, torch, os
from torch_geometric.data import extract_tar, download_url
from torch_pdb.datasets import TorchPDBDataset

class PDBBindRefined(TorchPDBDataset):
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
        return glob.glob(f'{self.root}/raw/files/*/*_protein.pdb')

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_refined.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_refined.tar.gz', f'{self.root}/raw')
        os.rename(f'{self.root}/raw/refined-set', f'{self.root}/raw/files')

    def add_protein_attributes(self, protein):
        pocket = self.pdb2df(f'{self.root}/raw/files/{protein["ID"]}/{protein["ID"]}_pocket.pdb')
        is_site = torch.zeros((len(pocket),))
        is_site[(
            torch.tensor(pocket['residue_number'].tolist()).unsqueeze(1) == protein['residue_index']
        ).sum(dim=1).nonzero()] = 1.
        protein['binding_site'] = is_site
        return protein
    def describe(self):
        desc = super().describe()
        desc['property'] = "Small Mol. Binding Site (`binding_site`)"
        desc['values'] = 2
        desc['type'] = 'Binary'
