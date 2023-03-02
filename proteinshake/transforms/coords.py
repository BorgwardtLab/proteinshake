from scipy.ndimage import rotate
import numpy as np

from proteinshake.transforms import Transform


def _get_coords_array(protein, resolution='residue'):
    """ Get a numpy array of the protein coordinates


    Arguments
    ---------
    protein: dict
        A protein dictionary.

    resolution: str
        The resolution at which to return the coordinates ('residue', or 'atom')

    Returns
    --------
    np.array
        Numpy array of shape (n_(residues or atoms),  3)

    """

    return np.array([protein['residue']['x'],
                     protein['residue']['y'],
                     protein['residue']['z']
                     ]
                   ).T

def _set_coords(protein, coord_array, resoluiton='residue', in_place=True):
    """ Given an Nx3 array of coordinates, set them to the
    protein dictionary coord key.

    Arguments
    ----------
    protein: dict
        Proteinshake protein dict
    resolution: str
        Which resolution to use ('residue' or 'atom')
    coord_array: np.array
        Nx3 numpy coordinate array

    """

    assert len(protein([resolution]['x']) == len(coord_array),\
            "Mismatch of coordinate array shape and protein, check choice of resolution"

    protein[resolution]['x'] = list(coord_array[:,0])
    protein[resolution]['y'] = list(coord_array[:,1])
    protein[resolution]['z'] = list(coord_array[:,2])

class CenterTransform(Transform):
    """ Center the coordinates of a protein at atom and residue level.
    We use the Ca to compute the center for all the atoms.
        """
    def __call__(self, protein):
        coords_res = _get_coords_array(protein, resolution='residue')
        coords_atom = _get_coords_array(protein, resolution='atom')

        center = np.mean(coords_res, axis=0)
        new_coords_res = coords_res - center
        new_coords_atom = coords_atom - center

        _set_coords(protein, new_coords_res, resolution='residue')
        _set_coords(protein, new_coords_atom, resolution='atom')

        return protein


class RandomRotateTransform(Transform):
    """ Apply a random rotation to the coordinate arrays of a protein"""

    def __init__(self, seed=42):
        self.seed = seed
        np.random.seed(self.seed)
        pass

    def __call__(self, protein):
        coords_res = _get_coords_array(protein, resolution='residue')
        coords_atom = _get_coords_array(protein, resolution='atom')

        # choose an angle and axis along which to rotate
        angle = np.random.randint(-180, 180, size=1)[0]
        axes_choice = [(0, 1), (0, 2), (1, 2)]
        axes = axes_choice[np.random.randint(0, len(axes_choices), size=1)[0]]

        coords_res_rot = rotate(coords_res, angle=angle, axes=axes)
        coords_atom_rot = rotate(coords_res, angle=angle, axes=axes)

        _set_coords(protein, coords_res_rot, resolution='residue')
        _set_coords(protein, coord_atom_rot, resolution='atom')

        return protein
