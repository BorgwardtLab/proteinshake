'''
Tests loading a dataset into different frameworks and the specific data loaders.
'''

import unittest, tempfile
from proteinshake.datasets import AlphaFoldDataset

class TestFrameworks(unittest.TestCase):

    def setUp(self):
        organism = 'methanocaldococcus jannaschii'
        self.tmpdir = tempfile.TemporaryDirectory()
        self.ds = AlphaFoldDataset(root=self.tmpdir, organism=organism)

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_graph_pyg(self):
        from torch_geometric.loader import DataLoader
        graphs = self.ds.to_graph(k=5).pyg()
        loader = DataLoader(graphs)
        x = next(iter(loader))

    def test_voxel_torch(self):
        from torch.utils.data import DataLoader
        voxels = self.ds.to_voxel().torch()
        loader = DataLoader(voxels)
        x = next(iter(loader))

    def test_point_torch(self):
        from torch.utils.data import DataLoader
        points = self.ds.to_point().torch()
        loader = DataLoader(points)
        x = next(iter(loader))


if __name__ == '__main__':
    unittest.main()
