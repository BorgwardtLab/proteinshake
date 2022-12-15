import os.path as osp
import tempfile
import shutil
import subprocess

import pandas as pd

""" Wrappers for external programs. """

def tmalign_wrapper(pdb1, pdb2):
    """Compute TM score with TMalign between two PDB structures.
    Parameters
    ----------
    pdb1: str
        Path to PDB.
    arg2 : str
        Path to PDB.
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
        out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return -1.
    return float(TM1), float(TM2), float(RMSD)


def cdhit_wrapper(sequences, sim_thresh=0.6):
    """ Cluster sequences using CD-hit

    Parameters
    -----------
    sequences: list
        List of protein sequences to cluster.

    Returns
    --------
    representatives: list
        List of sequence indices to preserve as representatives.
    """

    assert shutil.which('cd-hit') is not None,\
    "CD-HIT installation not found. Go here https://github.com/weizhongli/cdhit to install"

    with tempfile.TemporaryDirectory() as tmpdir:
        in_file = osp.join(tmpdir, 'in.fasta')
        out_file = osp.join(tmpdir, 'out.fasta')
        with open(in_file, "w") as inp:
            for i, s in enumerate(sequences):
                inp.write(f"> {i} \n")
                inp.write(s + "\n")

        try:
            cmd = ['cd-hit',
                   '-c',
                   str(sim_thresh),
                   '-i',
                   in_file,
                   '-o',
                   out_file
                  ]

            subprocess.run(cmd,
                           stdout=subprocess.DEVNULL,
                           stderr=subprocess.STDOUT
                          )
        except Exception as e:
            print(e)
            return -1
        else:
            clusters = [0] * len(sequences)
            with open(out_file + ".clstr", "r") as out:
                inds = []
                for line in out:
                    if line.startswith(">"):
                        clust_id = int(line.split()[1])
                        continue
                    ind = int(line.split(">")[1].split('.')[0])
                    clusters[ind] = clust_id
            return clusters

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
