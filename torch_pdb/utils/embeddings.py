"""
Some basic functions for embedding protein sequences. Supply to a representation class.

All embeddings take a sequence string as an input and return the embedding as a numpy array of shape n x d.
"""

import numpy as np

alphabet = 'ARNDCEQGHILKMFPSTWYV'

def one_hot(sequence):
    """ Compute the one-hot encoding of a protein sequence.

    Parameters
    ----------
    sequence: str
        The protein sequence

    Returns
    -------
    ndarray
        The embedded sequence.
    """

    return np.stack([np.eye(len(alphabet))[alphabet.index(aa)] for aa in sequence])

# AA Index
# ...

# ESM
# ...
