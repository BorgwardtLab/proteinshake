'''
Tests downloading the precomputed datasets and loading the data.
'''

import unittest, tempfile
from proteinshake.datasets import *
from proteinshake.datasets.alphafold import AF_DATASET_NAMES

class TestDatasets(unittest.TestCase):

    def test_pli(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandInterfaceDataset(root=tmp, verbosity=0).download_precomputed()

    def test_ppi(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinProteinInterfaceDataset(root=tmp, verbosity=0).download_precomputed()

    def test_tm(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = TMAlignDataset(root=tmp, verbosity=0).download_precomputed()

    def test_go(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = GeneOntologyDataset(root=tmp, verbosity=0).download_precomputed()

    def test_ec(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = EnzymeCommissionDataset(root=tmp, verbosity=0).download_precomputed()

    def test_pfam(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinFamilyDataset(root=tmp, verbosity=0).download_precomputed()

    def test_rcsb(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = RCSBDataset(root=tmp, verbosity=0).download_precomputed()

    def test_scop(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = SCOPDataset(root=tmp, verbosity=0).download_precomputed()

    def test_af(self):
        organism = 'methanocaldococcus jannaschii'
        with tempfile.TemporaryDirectory() as tmp:
            ds = AlphaFoldDataset(root=tmp, organism=organism, verbosity=0).download_precomputed()

    def test_dude(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandDecoysDataset(root=tmp, verbosity=0).download_precomputed()

if __name__ == '__main__':
    unittest.main()
