import shutil
import os
import subprocess
from biopandas.pdb import PandasPdb

def tmalign_wrapper(pdb1, pdb2, return_superposition=False):
    """Compute TM score with TMalign between two PDB structures.
    Parameters
    ----------
    pdb1: str
        Path to PDB.
    arg2 : str
        Path to PDB.
    return_superposition: bool
        If True, returns a protein dataframe with superposed structures.
    Returns
    -------
    float
        TM score from `pdb1` to `pdb2`
    float
        TM score from `pdb2` to `pdb1`
    float
        RMSD between structures
    """
    assert shutil.which('TMalign') is not None,\
           "No TMalign installation found. Go here to install : https://zhanggroup.org/TM-align/TMalign.cpp"
    try:
        if return_superposition:
            with tempfile.TemporaryDirectory() as tmpdir:
                out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2, '-o', os.path.join(tmpdir, 'superposition.pdb')], stdout=subprocess.PIPE).stdout.decode()
                superposition = PandasPdb().read_pdb(os.path.join(tmpdir, 'superposition.pdb')) 
        else:
            out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return -1.
    if return_superposition:
        return float(TM1), float(TM2), float(RMSD), superposition
    else:
        return float(TM1), float(TM2), float(RMSD)