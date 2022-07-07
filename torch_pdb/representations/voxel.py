from tqdm import tqdm

class VoxelDataset():

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
