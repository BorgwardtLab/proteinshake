import unittest, tempfile

from proteinshake.datasets import ProteinLigandDecoysDataset as dset
from proteinshake.transforms import *


class TestTransforms(unittest.TestCase):

    def test_identity(self):
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_voxel(transform=IdentityTransform())

    def test_center(self):
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_voxel(transform=CenterTransform())

    def test_rotate(self):
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_voxel(transform=RandomRotateTransform())

    def test_compose(self):
        t = Compose([CenterTransform(), RandomRotateTransform()])
        with tempfile.TemporaryDirectory() as tmp:
            da = dset(root=tmp).to_voxel(transform=t)
