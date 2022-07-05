import torch

alphabet = 'ARNDCEQGHILKMFPSTWYV'

def one_hot(sequence):
    return torch.stack([torch.eye(len(alphabet))[alphabet.index(aa)] for aa in sequence])

# AA Index
# ...

# ESM
# ...
