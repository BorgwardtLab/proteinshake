import os
import subprocess

from tqdm import tqdm
import numpy as np
import open3d as o3d

from proteinshake.utils import protein_to_pdb

class Surface():
    """ Surface representation of a protein.

    Converts a protein object to a surface using open3d.

    Parameters
    ----------
    protein: dict
        A protein object.
    construction: str
        Whether to use knn or eps construction.

    """

    def __init__(self, protein):
        resolution = 'atom' if 'atom' in protein else 'residue'

        self.data = self._compute_surface(protein)

    def _compute_surface(self, protein):
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
            protein_to_pdb(protein, pdb_path)
            cmd = ['dms', pdb_path, '-n', '-o', 'surf.txt']
            pass

    def _parse_dms(self, path):
        """
        LYS 215A    N-100.374   29.066   -7.068 A
        LYS 215A    N-101.056   27.395   -7.760 SR0  0.185 -0.756 -0.000 -0.655
        LYS 215A    N-101.148   27.825   -7.760 SR0  0.185 -0.690 -0.307 -0.655
        LYS 215A    N-100.760   28.179   -8.265 SR0  0.185 -0.437 -0.318 -0.841
        """
        pass

class SurfaceDataset():
    """ Graph representation of a protein structure dataset.

    Parameters
    ----------
    proteins: generator
        A generator of protein objects from a Dataset.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    resolution: str, default 'residue'
        Resolution of the proteins to use in the graph representation. Can be 'atom' or 'residue'.
    method: str
        Surface reconstruction algorithms (see http://www.open3d.org/docs/latest/tutorial/geometry/surface_reconstruction.html).
        Can be 'ball', 'alpha', 'poisson', Default is 'ball'.
    """

    def __init__(self, proteins, size, path, method='ball', resolution='residue'):
        self.path = f'{path}/processed/surface/{resolution}_{method}'
        self.surfaces = (Surface(protein) for protein in proteins)
        self.size = size
        os.makedirs(os.path.dirname(self.path), exist_ok=True)

    def pyg(self, *args, **kwargs):
        from proteinshake.frameworks.pyg import PygGraphDataset
        return PygGraphDataset(self.surfaces, self.size, self.path+'.pyg', *args, **kwargs)

    def dgl(self, *args, **kwargs):
        from proteinshake.frameworks.dgl import DGLGraphDataset
        return DGLGraphDataset(self.surfaces, self.size, self.path+'.dgl', *args, **kwargs)

    def nx(self, *args, **kwargs):
        from proteinshake.frameworks.nx import NetworkxGraphDataset
        return NetworkxGraphDataset(self.surfaces, self.size, self.path+'.nx', *args, **kwargs)

if __name__ == "__main__":
    from proteinshake.datasets import TMAlignDataset
    import tempfile
    with tempfile.TemporaryDirectory() as tf:
        da = TMAlignDataset(root=tf)
        da_surf = da.to_surface()
        list(da_surf.surfaces)
