import os
import shutil
import subprocess
import tempfile

import pandas as pd

from proteinshake.transforms import ShakeTransform
from proteinshake.utils import protein_to_pdb

class SurfaceTransform(ShakeTransform):
    def __init__(self, d=0.2):
        super().__init__()
        self.d = d

    def __call__(self, protein):
        surf = self._compute_surface(protein)
        mode = 'atom' if 'atom' in protein.keys() else 'residue'
        new_prot = dict()
        new_prot['protein'] = protein['protein']
        new_prot['protein']['residue'] = self._surf_filter(protein, surf, resolution='residue')
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
        print(surface_df.head())
        new_data = {
            'residue_number': surface_df['residue_index'].tolist(),
            'residue_type': surface_df['residue_name'].tolist(),
            'x': surface_df['x'].tolist(),
            'y': surface_df['y'].tolist(),
            'z': surface_df['z'].tolist(),
            'x_norm': surface_df['x_norm'].tolist(),
            'y_norm': surface_df['y_norm'].tolist(),
            'z_norm': surface_df['z_norm'].tolist(),
            }

        if resolution == 'atom':
            new_data['atom_number'] = list(range(len(surface_df))),
            new_data['atom_type'] = surface_df['atom_name'].tolist(),
        return new_data

    def _compute_surface(self, protein, d=0.2):
        """ Call DMS to compute a surface for the PDB.

        Usage: dms input_file [-a] [-d density] [-g file] [-i file] [-n] [-w radius] [-v] -o file
        -a	use all atoms, not just amino acids
        -d	change density of points
        -g	send messages to file
        -i	calculate only surface for specified atoms
        -n	calculate normals for surface points
        -w	change probe radius
        -v	verbose
        -o	specify output file name (required)

        See: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html#ref
        """
        with tempfile.TemporaryDirectory() as tf:
            pdb_path = os.path.join(tf, "in.pdb")
            dest = os.path.join(tf, "out.surf")
            protein_to_pdb(protein, pdb_path)
            assert shutil.which('dms') is not None, "DMS executable not in PATH go here to install https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/dms1.html#ref."
            cmd = ['dms', pdb_path, '-n', '-d', str(d), '-o', dest]
            subprocess.run(cmd,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT
                           )
            return self._parse_dms(dest)

    def _parse_dms(self, path):
        """ Extract surface points and normal vectors for each
        point.

        Parameters
        ------------
        path:
            Path to DMS output file.

        Returns
        --------
        pd.DataFrame
            DataFrame with one row for each surface atom.
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
        df = pd.read_csv(path,
                         delim_whitespace=True,
                         header=None,
                         names=names
                         )
        df = df.dropna(axis=0)
        return df

if __name__ == "__main__":
    from proteinshake.datasets import TMAlignDataset
    from proteinshake.transforms import SurfaceTransform
    import tempfile
    with tempfile.TemporaryDirectory() as tf:
        da = TMAlignDataset(root=tf)
        da_surf = da.to_point(
                              resolution='atom',
                              transform=SurfaceTransform()
                              ).torch()
