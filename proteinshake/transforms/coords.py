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

def _set_coords(protein, coord_array, resolution='residue', in_place=True):
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

    protein[resolution]['x'] = list(coord_array[:,0])
    protein[resolution]['y'] = list(coord_array[:,1])
    protein[resolution]['z'] = list(coord_array[:,2])

class CenterTransform(Transform):
    """ Center the coordinates of a protein at atom and residue level.
    We use the Ca to compute the center for all the atoms.
        """
    def __init__(self, resolution='residue'):
        self.resolution = resolution
        super().__init__()

    def __call__(self, protein):
        coords = _get_coords_array(protein, resolution=self.resolution)

        center = np.mean(coords, axis=0)
        new_coords = coords - center

        _set_coords(protein, new_coords, resolution=self.resolution)

        return protein

class RandomRotateTransform(Transform):
    """ Apply a random rotation to the coordinate arrays of a protein"""

    def __init__(self, resolution='residue', seed=42):
        self.seed = seed
        self.resolution = resolution
        np.random.seed(self.seed)
        super().__init__()
        pass

    def __call__(self, protein):
        coords = _get_coords_array(protein, resolution=self.resolution)

        # choose an angle and axis along which to rotate
        angle = np.random.randint(-180, 180, size=1)[0]
        axes_choice = [(0, 1), (0, 2), (1, 2)]
        axes = axes_choice[np.random.randint(0, len(axes_choice), size=1)[0]]
        print(axes)

        coords_rot = rotate(coords, angle=angle, axes=axes)

        _set_coords(protein, coords_rot, resolution=self.resolution)

        return protein
