'''
Tests loading a dataset into different frameworks and the specific data loaders.
'''

import unittest, tempfile
from proteinshake.datasets import ProteinLigandDecoysDataset

class TestFrameworks(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.ds = ProteinLigandDecoysDataset(root=self.tmpdir)

    @classmethod
    def tearDownClass(self):
        self.tmpdir.cleanup()

    def test_graph_pyg(self):
        from torch_geometric.loader import DataLoader
        graphs = self.ds.to_graph(k=5).pyg()
        loader = DataLoader(graphs)
        x = next(iter(loader))

    def test_graph_dgl(self):
        from dgl.dataloading import GraphDataLoader as DataLoader
        graphs = self.ds.to_graph(k=5).dgl()
        loader = DataLoader(graphs)
        x = next(iter(loader))

    def test_graph_nx(self):
        graphs = self.ds.to_graph(k=5).nx()
        x = graphs[0]

    def test_voxel_torch(self):
        from torch.utils.data import DataLoader
        voxels = self.ds.to_voxel().torch()
        loader = DataLoader(voxels)
        x = next(iter(loader))

    def test_voxel_tf(self):
        import tensorflow as tf
        voxels = self.ds.to_voxel().tf()
        def generator():
            for data, protein_dict in voxels:
                yield data
        loader = tf.data.Dataset.from_generator(generator, output_types=(tf.float32))
        x = next(iter(loader))

    def test_voxel_np(self):
        import numpy as np
        points = self.ds.to_point().np()
        def generator():
            for data, protein_dict in points:
                yield data
        loader = np.fromiter(generator(), object)
        x = next(iter(loader))

    def test_point_torch(self):
        from torch.utils.data import DataLoader
        points = self.ds.to_point().torch()
        loader = DataLoader(points)
        x = next(iter(loader))

    def test_point_tf(self):
        import tensorflow as tf
        points = self.ds.to_point().tf()
        def generator():
            for data, protein_dict in points:
                yield data
        loader = tf.data.Dataset.from_generator(generator, output_types=(tf.float32))
        x = next(iter(loader))

    def test_point_np(self):
        import numpy as np
        points = self.ds.to_point().np()
        def generator():
            for data, protein_dict in points:
                yield data
        loader = np.fromiter(generator(), object)
        x = next(iter(loader))


if __name__ == '__main__':
    unittest.main()
