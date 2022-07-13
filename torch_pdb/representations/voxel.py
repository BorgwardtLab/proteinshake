from tqdm import tqdm

class VoxelDataset():
    """ Voxel representation of a protein structure dataset.

    Voxelizes the protein dataset.

    Parameters
    ----------
    embedding: Union[function, list]
        A function or list of functions for embedding the protein sequence to node attributes.

    """

    def __init__(self, root, proteins):
        self.root = root
        self.name = f'x'
        self.proteins = self.convert(proteins)

    def protein2voxel(self, protein):
        pass

    def voxel2torch(self, graph, info={}):
        pass

    def convert(self, proteins):
        return [self.protein2voxel(p) for p in tqdm(proteins, desc='Voxelizing Proteins')]

    def torch(self):
        raise NotImplementedError

    def tf(self):
        raise NotImplementedError
