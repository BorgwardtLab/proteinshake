import shutil
import os
import subprocess
import tempfile
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
import numpy as np

def global_distance_test(superposition):
    df = superposition.df['ATOM']
    df = df[df['atom_name'] == 'CA']
    A = df[df['chain_id'] == 'A']
    B = df[df['chain_id'] == 'B']
    coordsA = np.array(list(zip(A['x_coord'], A['y_coord'], A['z_coord'])))
    coordsB = np.array(list(zip(B['x_coord'], B['y_coord'], B['z_coord'])))
    dist = np.sqrt(((coordsA-coordsB)**2).sum(1))
    thresholds = [1,2,4,8]
    return np.mean([np.mean(dist<=t) for t in thresholds])

def tmalign_wrapper(pdb1, pdb2):
    """Compute TM score with TMalign between two PDB structures.
    Parameters
    ----------
    pdb1: str
        Path to PDB.
    pdb2 : str
        Path to PDB.
    return_superposition: bool
        If True, returns a protein dataframe with superposed structures.
    Returns
    -------
    dict
        Metric values TM1/TM2 (TM-Scores normalized to pdb1 or pdb2), RMSD, GDT
    """
    assert shutil.which('TMalign') is not None,\
           "No TMalign installation found. Go here to install : https://zhanggroup.org/TM-align/TMalign.cpp"
    with tempfile.TemporaryDirectory() as tmpdir:
        out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2, '-o', os.path.join(tmpdir, 'superposition.pdb')], stdout=subprocess.PIPE).stdout.decode()
        superposition = PandasPdb().read_pdb(os.path.join(tmpdir, 'superposition.pdb')) 
    path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    GDT = global_distance_test(superposition)
    return {
        'TM1': float(TM1),
        'TM2': float(TM2),
        'RMSD': float(RMSD),
        'GDT': GDT,
    }