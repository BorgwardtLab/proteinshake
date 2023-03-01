"""
Abstract class for transforming a protein.
"""

class Transform:
    def __init__(self):
        pass
    def __call__(self, protein):
        """ Accepts a protein dict and returns a new protein dict
        """
        raise NotImplementedError

class IdentityTransform(Transform):
    def __call__(self, protein):
        return protein
