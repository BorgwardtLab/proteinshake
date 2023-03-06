# -*- coding: utf-8 -*-
import glob
import os

from biopandas.pdb import PandasPdb

from proteinshake.datasets import Dataset
from proteinshake.utils import download_url

AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'CYX': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HID': 'H', 'HIE': 'HIS', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

class ProteinLigandDecoyDataset(Dataset):
    """ Proteins (targets) from `DUDE-Z <https://pubs.acs.org/doi/10.1021/acs.jcim.0c00598>`_ with a list of decoys and active molecules for each.
    Each molecule is encoded as a SMILES string, meant to be used in a virtual screen setting.
    In this setting a model is given a protein and a ligand and outputs a score reflecting the likelihood
    that the given molecule is a binder. Then, this score is used to sort the union of all the ligands and
    decoys. A good model places true ligands at the top of this list. This is known as enrichment factor
    analysis.


    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.

    .. code-block:: python

        >>> from proteinshake.datasets import ProteinLigandDecoyDataset
        >>> dataset = ProteinLigandDecoyDataset()
        >>> protein = next(dataset.proteins())
        >>> protein['ligands'][:3]
        >>> protein['decoys'][:3]

    """

    def pdb2df(self, path):
        """ Parses a single PDB file to a DataFrame (with biopandas).
        Also deals with multiple structure models in a PDB (e.g. from NMR) by only selecting the first model.
        This modification allows for protonated histidine residues which are common in this dataset.
        We just map them to regular HIS for now.

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
        df = PandasPdb().read_pdb_from_list(filtered_lines).df['ATOM']
        print(df['residue_name'].unique())
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



    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/*.pdb')[:self.limit]

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):

        targets  = [
                    'AA2AR',
                    'ABL1',
                    'ACES',
                    'ADA',
                    'ADRB2',
                    'AMPC',
                    'ANDR',
                    'CSF1R',
                    'CXCR4',
                    'DEF',
                    'DRD4',
                    'EGFR',
                    'FA7',
                    'FA10',
                    'FABP4',
                    'FGFR1',
                    'FKB1A',
                    'GLCM',
                    'HDAC8',
                    'HIVPR',
                    'HMDH',
                    'HS90A',
                    'ITAL',
                    'KITH',
                    'KIT',
                    'LCK',
                    'MAPK2',
                    'MK01',
                    'MT1',
                    'NRAM',
                    'PARP1',
                    'PLK1',
                    'PPARA',
                    'PTN1',
                    'PUR2',
                    'RENI',
                    'ROCK1',
                    'SRC',
                    'THRB',
                    'TRY1',
                    'TRYB1',
                    'UROK',
                    'XIAP',
                ]


        for target_id in targets[:3]:
            # grab receptor
            print(target_id)
            download_url(f"https://dudez.docking.org/DOCKING_GRIDS_AND_POSES/{target_id}/rec.crg.pdb",
                         f"{self.root}/raw/"
                        )
            os.rename(f'{self.root}/raw/rec.crg.pdb', f'{self.root}/raw/{target_id}.pdb')
            # grab ligands
            download_url(f"https://dudez.docking.org/property_matched/{target_id}_new_DUDE_1/ligands.smi",
                        f"{self.root}/raw/"
                        )
            os.rename(f'{self.root}/raw/ligands.smi', f'{self.root}/raw/ligands_{target_id}.smi')

            # grab decoys
            download_url(f"https://dudez.docking.org/property_matched/{target_id}_new_DUDE_1/decoys.smi",
                         f"{self.root}/raw/"
                        )
            os.rename(f'{self.root}/raw/decoys.smi', f'{self.root}/raw/decoys_{target_id}.smi')


    def add_protein_attributes(self, protein):

        target = protein['protein']['ID']
        with open(f"{self.root}/raw/decoys_{target}.smi", "r") as ligands:
            protein['protein']['ligands_smiles'] = [line.split()[0] for line in ligands.readlines()]
            protein['protein']['ligands_ids'] = [line.split()[1] for line in ligands.readlines()]

        with open(f"{self.root}/raw/decoys_{target}.smi", "r") as decoys:
            protein['protein']['decoys_smiles'] = [line.split()[0] for line in decoys.readlines()]
            protein['protein']['decoys_ids'] = [line.split()[1] for line in decoys.readlines()]

        return protein

if __name__ == "__main__":
    da = ProteinLigandDecoyDataset(use_precomputed=False)
    protein = next(da.proteins())
    # print(protein)

