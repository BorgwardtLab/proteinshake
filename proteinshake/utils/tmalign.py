import shutil
import os
import subprocess
import tempfile
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
import numpy as np

def global_distance_test(superposition, alignmentA, alignmentB):
    df = superposition.df['ATOM']
    A = df[df['chain_id'] == 'A']
    B = df[df['chain_id'] == 'B']
    coordsA = np.array(list(zip(A['x_coord'], A['y_coord'], A['z_coord'])))[alignmentA]
    coordsB = np.array(list(zip(B['x_coord'], B['y_coord'], B['z_coord'])))[alignmentB]
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
        lines = subprocess.run(['TMalign','-outfmt','-1', pdb1, pdb2, '-o', f'{tmpdir}/superposition'], stdout=subprocess.PIPE).stdout.decode().split('\n')
        TM1 = lines[7].split()[1]
        TM2 = lines[8].split()[1]
        RMSD = lines[6].split()[4][:-1]
        seq1, ali, seq2 = lines[12], lines[13], lines[14]
        i, j, alignmentA, alignmentB = 0, 0, [], []
        for s1,a,s2 in zip(seq1,ali,seq2):
            if a != ' ': alignmentA.append(i)
            if a != ' ': alignmentB.append(j)
            if s1 != '-': i += 1
            if s2 != '-': j += 1
        os.rename(f'{tmpdir}/superposition_all', f'{tmpdir}/superposition_all.pdb')
        superposition = PandasPdb().read_pdb(f'{tmpdir}/superposition_all.pdb')
        GDT = global_distance_test(superposition, alignmentA, alignmentB)
    return {
        'TM1': float(TM1),
        'TM2': float(TM2),
        'RMSD': float(RMSD),
        'GDT': GDT,
    }