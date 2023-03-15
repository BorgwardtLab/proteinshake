import unittest, tempfile

from proteinshake.datasets import ProteinLigandDecoysDataset

from proteinshake.transforms import *


class TestTransforms(unittest.TestCase):

    @classmethod
    def setUpClass(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        self.ds = ProteinLigandDecoysDataset(root=self.tmpdir)

    @classmethod
    def tearDownClass(self):
        self.tmpdir.cleanup()

    def test_identity(self):
        self.ds.to_voxel(transform=IdentityTransform())

    def test_center(self):
        self.ds.to_voxel(transform=CenterTransform())

    def test_rotate(self):
        self.ds.to_voxel(transform=RandomRotateTransform())

    def test_compose(self):
        self.ds.to_voxel(transform=Compose([CenterTransform(), RandomRotateTransform()]))
