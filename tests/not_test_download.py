'''
Tests all downloads with 'use_precomputed=False'. The number of downloaded files and the number of parsed files is patched to a small number where possible. However, most datasets require downloading one large file, which takes time. Hence removed from GitHub testing CI (by not naming it according to pytest convention).
'''

import unittest, tempfile
from unittest import mock
from proteinshake.datasets import (RCSBDataset,
                                   GeneOntologyDataset,
                                   EnzymeCommissionDataset,
                                   PfamDataset,
                                   ProteinProteinInterfaceDataset,
                                   ProteinLigandInterfaceDataset,
                                   TMAlignDataset,
                                   AlphaFoldDataset,
                                   SCOPDataset
                                   )
from proteinshake.datasets.alphafold import AF_DATASET_NAMES

class TestDownload(unittest.TestCase):

    @mock.patch('proteinshake.datasets.ProteinLigandInterfaceDataset.download_limit', return_value=5)
    def test_protein_ligand(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandInterfaceDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.ProteinProteinInterfaceDataset.download_limit', return_value=5)
    def test_protein_protein(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinProteinInterfaceDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.RCSBDataset.download_limit', return_value=5)
    def test_rcsb(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = RCSBDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.GeneOntologyDataset.download_limit', return_value=5)
    def test_go(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = GeneOntologyDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.EnzymeCommissionDataset.download_limit', return_value=5)
    def test_go(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = EnzymeCommissionDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.PfamDataset.download_limit', return_value=5)
    def test_pfam(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = PfamDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.TMAlignDataset.download_limit', return_value=5)
    def test_pfam(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = TMAlignDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.AlphaFoldDataset.download_limit', return_value=5)
    def test_pfam(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            organism = 'methanocaldococcus jannaschii'
            ds = AlphaFoldDataset(root=tmp, organism=organism, use_precomputed=False)

    @mock.patch('proteinshake.datasets.SCOPDataset.download_limit', return_value=5)
    def test_scop(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = SCOPDataset(root=tmp, use_precomputed=False)

if __name__ == '__main__':
    unittest.main()
