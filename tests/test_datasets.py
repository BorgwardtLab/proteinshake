'''
Tests downloading the precomputed datasets and loading the data.
'''

import unittest, tempfile
from proteinshake.datasets import *
from proteinshake.datasets.alphafold import AF_DATASET_NAMES

class TestDatasets(unittest.TestCase):

    def test_pli(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandInterfaceDataset(root=tmp).download_precomputed()

    def test_ppi(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinProteinInterfaceDataset(root=tmp).download_precomputed()

    def test_tm(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = TMAlignDataset(root=tmp).download_precomputed()

    def test_go(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = GeneOntologyDataset(root=tmp).download_precomputed()

    def test_ec(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = EnzymeCommissionDataset(root=tmp).download_precomputed()

    def test_pfam(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinFamilyDataset(root=tmp).download_precomputed()

    def test_rcsb(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = RCSBDataset(root=tmp).download_precomputed()

    def test_af(self):
        organism = 'methanocaldococcus jannaschii'
        with tempfile.TemporaryDirectory() as tmp:
            ds = AlphaFoldDataset(root=tmp, organism=organism).download_precomputed()

    def test_dude(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandDecoysDataset(root=tmp).download_precomputed()

if __name__ == '__main__':
    unittest.main()
