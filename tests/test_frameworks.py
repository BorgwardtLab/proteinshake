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

    def test_pyg(self):
        from torch_geometric.loader import DataLoader
        graphs = self.ds.to_graph(k=5).pyg()
        loader = DataLoader(graphs)
        x = next(iter(loader))


if __name__ == '__main__':
    unittest.main()
