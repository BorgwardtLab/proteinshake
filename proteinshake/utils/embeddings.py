"""
Some basic functions for embedding protein sequences. Supply to a representation class.

All embeddings take a sequence string as an input and return the embedding as a numpy array of shape n x d.
"""

import numpy as np

residue_alphabet = 'ARNDCEQGHILKMFPSTWYV'
atom_alphabet = 'NCOSH'

def onehot(sequence, resolution='residue'):
    """ Compute the one-hot encoding of a protein sequence.

    Parameters
    ----------
    sequence: str
        The protein sequence.
    resolution: str, default 'resolution'
        Resolution of the protein. 'residue' or 'atom'.

    Returns
    -------
    ndarray
        The embedded sequence.
    """
    if resolution == 'residue':
        return np.stack([np.eye(len(residue_alphabet))[residue_alphabet.index(aa)] for aa in sequence])
    else:
        return np.stack([np.eye(len(atom_alphabet))[atom_alphabet.index(aa[0])] for aa in sequence])



def tokenize(sequence, resolution='residue'):
    """ Tokenizes the sequence.

    Parameters
    ----------
    sequence: str
        The protein sequence.
    resolution: str, default 'resolution'
        Resolution of the protein. 'residue' or 'atom'.

    Returns
    -------
    ndarray
        The embedded sequence.
    """
    if resolution == 'residue':
        return np.array([residue_alphabet.index(aa) for aa in sequence])
    else:
        return np.array([atom_alphabet.index(aa[0]) for aa in sequence])

# from: https://gist.github.com/foowaa/5b20aebd1dff19ee024b6c72e14347bb
def sinusoid_encoding_table(n_position, d_hid, padding_idx=None):
    """ Helper function to build a sinusoidal lookup table for positional encodings.
    """
    def cal_angle(position, hid_idx):
        return position / np.power(10000, 2 * (hid_idx // 2) / d_hid)
    def get_posi_angle_vec(position):
        return [cal_angle(position, hid_j) for hid_j in range(d_hid)]
    sinusoid_table = np.array([get_posi_angle_vec(pos_i) for pos_i in range(n_position)])
    sinusoid_table[:, 0::2] = np.sin(sinusoid_table[:, 0::2])
    sinusoid_table[:, 1::2] = np.cos(sinusoid_table[:, 1::2])
    if padding_idx is not None:
        sinusoid_table[padding_idx] = 0.
    return sinusoid_table

def positional_encoding(sequence, dim=128):
    """ Sinusoidal encoding of sequence position.

    Parameters
    ----------
    sequence: str
        The protein sequence

    Returns
    -------
    ndarray
        The embedded sequence.
    """
    n = len(sequence)
    table = sinusoid_encoding_table(n, dim)
    return table

def compose_embeddings(embeddings):
    """ Composes multiple embeddings into one by concatenating the results.

    Parameters
    ----------
    embeddings: list
        A list of embeddings

    Returns
    -------
    function
        A substitute embedding function.
    """

    return lambda sequence: np.concatenate([e(sequence) for e in embeddings], axis=-1)
