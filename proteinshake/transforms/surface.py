import os
import shutil
import subprocess
import tempfile

import pandas as pd

from proteinshake.transforms import ShakeTransform
from proteinshake.utils import protein_to_pdb
from proteinshake.utils import dms_wrapper

AA_THREE_TO_ONE = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

class SurfaceTransform(ShakeTransform):
    def __init__(self, d=0.2):
        super().__init__()
        self.d = d

    def __call__(self, protein):
        surf = dms_wrapper(protein)
        mode = 'atom' if 'atom' in protein.keys() else 'residue'
        new_prot = dict()
        new_prot['protein'] = protein['protein']
        new_prot['residue'] = self._surf_filter(protein, surf, resolution='residue')
        if mode == 'atom':
            new_prot['atom'] = self._surf_filter(protein, surf, resolution='atom')
        return new_prot

    def _surf_filter(self, protein, surface_df, resolution='residue'):

        """
        names = ['residue_name',
                 'residue_index',
                 'atom_name',
                 'x',
                 'y',
                 'z',
                 'point_type',
                 'area',
                 'x_norm',
                 'y_norm',
                 'z_norm'
                 ]
        """
        new_data = {
            'residue_number': surface_df['residue_index'].tolist(),
            'residue_type': [AA_THREE_TO_ONE[a] for a in surface_df['residue_name'].tolist()],
            'x': surface_df['x'].tolist(),
            'y': surface_df['y'].tolist(),
            'z': surface_df['z'].tolist(),
            'x_norm': surface_df['x_norm'].tolist(),
            'y_norm': surface_df['y_norm'].tolist(),
            'z_norm': surface_df['z_norm'].tolist(),
            }

        if resolution == 'atom':
            new_data['atom_number'] = list(range(len(surface_df)))
            new_data['atom_type'] = surface_df['atom_name'].tolist()
        
        return new_data
