# -*- coding: utf-8 -*-
import glob
import os
import os.path as osp

from rdkit import Chem
from rdkit import RDLogger
import numpy as np

from proteinshake.datasets import Dataset
from proteinshake.utils.pdbbind import parse_pdbbind_PL_index
from proteinshake.utils import extract_tar, download_url

RDLogger.DisableLog('rdApp.*') # disable warnings

class ProteinLigandInterfaceDataset(Dataset):
    """Proteins bound to small molecules with binding site and affinity information. Residues
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
        return glob.glob(f'{self.root}/raw/files/*/*_protein.pdb')[:self.download_limit()]

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_refined.tar.gz', f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_refined.tar.gz', f'{self.root}/raw')
        os.rename(f'{self.root}/raw/refined-set', f'{self.root}/raw/files')

    def add_protein_attributes(self, protein):
        pocket = self.pdb2df(f'{self.root}/raw/files/{protein["protein"]["ID"]}/{protein["protein"]["ID"]}_pocket.pdb')
        is_site = np.zeros((len(pocket),))
        is_site[(
            np.expand_dims(np.array(pocket['residue_number'].tolist()), axis=1) == protein['residue']['residue_number']
        ).sum(axis=1).nonzero()] = 1.
        protein['residue']['binding_site'] = is_site

        index_data = parse_pdbbind_PL_index(osp.join(self.root,
                                                    "raw",
                                                    "files",
                                                    "index",
                                                    f"INDEX_refined_set.{self.version}")
                                            )
        ligand = Chem.MolFromMolFile(osp.join(self.root,
                                              'raw',
                                              'files',
                                              protein['protein']['ID'],
                                              f'{protein["protein"]["ID"]}_ligand.sdf')
                                     )
        try:
            smiles = Chem.MolToSmiles(ligand)
        except:
            return None
        is_site = np.zeros((len(pocket),))
        is_site[(
            np.expand_dims(np.array(pocket['residue_number'].tolist()), axis=1) == protein['residue']['residue_number']
        ).sum(axis=1).nonzero()] = 1.
        protein['residue']['binding_site'] = is_site.tolist()
        bind_data = index_data[protein['protein']['ID']]
        protein['protein']['kd'] = bind_data['kd']['value']
        protein['protein']['resolution'] = bind_data['resolution']
        protein['protein']['year'] = bind_data['date']
        protein['protein']['ligand_id'] = bind_data['ligand_id']
        protein['protein']['ligand_smiles'] = smiles
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Small Mol. Binding Site (residue-level)"
        desc['values'] = 2
        desc['type'] = 'Binary'
        return desc
