"""
Abstract class for transforming a protein.
"""

class Transform:
    """ A callable object which accepts a protein dictionary and returns an updated version of it."""
    def __init__(self):
        pass

    def __call__(self, protein):
        """ Accepts a protein dict and returns a new protein dict

        Arguments
        ----------
        protein: dict
            A protein dictionary.

        """
        raise NotImplementedError

class IdentityTransform:
    """ Do nothing to the protein"""
    def __call__(self, protein):
        return protein

class Compose:

    def __init__(self, transforms):
        self.transforms = transforms

    def __call__(self, data):
        for transform in self.transforms:
            data = transform(data)
        return data

    def __repr__(self):
        args = [f'  {transform}' for transform in self.transforms]
        return '{}([\n{}\n])'.format(self.__class__.__name__, ',\n'.join(args))