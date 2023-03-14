import unittest, tempfile

from proteinshake.datasets import ProteinLigandDecoysDataset as dset


class TestRepresentations(unittest.TestCase):

    def test_graph(self):
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_voxel()

    def test_voxel(self):
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_graph(eps=9)

    def test_point(self):
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_point()
