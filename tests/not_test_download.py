'''
Tests all downloads with 'use_precomputed=False'. The number of downloaded files and the number of parsed files is patched to a small number where possible. However, most datasets require downloading one large file, which takes time. Hence removed from GitHub testing CI (by not naming it according to pytest convention).
'''

import unittest, tempfile
from unittest import mock
from torch_pdb.datasets import PDBBindRefined, PDBBindPPI,  TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset, AlphaFoldDataset
from torch_pdb.datasets.alphafold import AF_DATASET_NAMES

class TestDownload(unittest.TestCase):

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_pdbbind(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PDBBindRefined(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_pdbbind_ppi(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PDBBindPPI(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_tm(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = TMScoreBenchmark(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_go(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = GODataset(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_ec(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ECDataset(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_pfam(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PfamDataset(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_rcsb(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = RCSBDataset(root=tmp, use_precomputed=False)

    @mock.patch('torch_pdb.datasets.TorchPDBDataset.download_limit', return_value=5)
    def test_alphafold(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = AlphaFoldDataset(root=tmp, organism='methanocaldococcus jannaschii', use_precomputed=False)

if __name__ == '__main__':
    unittest.main()
