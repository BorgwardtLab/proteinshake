import os
import tempfile
import shutil
import subprocess
import pandas as pd

from proteinshake.utils import protein_to_pdb

def dms_wrapper(protein, d=0.2):
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
        return _parse_dms(dest)

def _parse_dms(path):
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