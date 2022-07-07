from tqdm import tqdm

class PointDataset():

    def __init__(self, root, proteins, embedding):
        self.root = root
        self.name = f'embedding_{embedding.__name__}'
        self.proteins = self.convert(proteins)

    def protein2point(self, protein):
        pass

    def point2torch(self, graph, info={}):
        pass

    def convert(self, proteins):
        return [self.protein2point(p) for p in tqdm(proteins, desc='Converting proteins to point clouds')]

    def torch(self):
        raise NotImplementedError

    def tf(self):
        raise NotImplementedError
