import torch

alphabet = 'ARNDCEQGHILKMFPSTWYV'

def one_hot(sequence):
    """ Compute the one-hot encoding of a protein sequence.

    Args:
        sequence:

    Returns:
    """

    return torch.stack([torch.eye(len(alphabet))[alphabet.index(aa)] for aa in sequence])

# AA Index
# ...

# ESM
# ...
