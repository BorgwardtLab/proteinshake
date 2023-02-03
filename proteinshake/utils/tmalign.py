import shutil
import os
import subprocess
import tempfile
from biopandas.pdb import PandasPdb
from .gdt import gdt

def tmalign_wrapper(pdb1, pdb2):
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
    with tempfile.TemporaryDirectory() as tmpdir:
        out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2, '-o', os.path.join(tmpdir, 'superposition.pdb')], stdout=subprocess.PIPE).stdout.decode()
        superposition = PandasPdb().read_pdb(os.path.join(tmpdir, 'superposition.pdb')) 
    path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    GDT = gdt(superposition)
    return {
        'TM1': float(TM1),
        'TM2': float(TM2),
        'RMSD': float(RMSD),
        'GDT': GDT,
    }