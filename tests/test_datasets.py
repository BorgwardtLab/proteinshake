'''
Tests downloading the precomputed datasets and loading the data.
'''

import unittest, tempfile
from torch_pdb.datasets import PDBBindRefined, PDBBindPPI, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset, AlphaFoldDataset
from torch_pdb.datasets.alphafold import AF_DATASET_NAMES

class TestDatasets(unittest.TestCase):

    def test_pdbbind(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PDBBindRefined(root=tmp)

    def test_pdbbind_ppi(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PDBBindPPI(root=tmp)

    def test_tm(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = TMScoreBenchmark(root=tmp)

    def test_go(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = GODataset(root=tmp)

    def test_ec(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ECDataset(root=tmp)

    def test_pfam(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PfamDataset(root=tmp)

    def test_rcsb(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = RCSBDataset(root=tmp)

    def test_alphafold(self):
        organism = 'methanocaldococcus jannaschii'
        with tempfile.TemporaryDirectory() as tmp:
            ds = AlphaFoldDataset(root=tmp, organism=organism)

if __name__ == '__main__':
    unittest.main()
