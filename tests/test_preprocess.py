'''
Tests all downloads with 'use_precomputed=False'. The number of downloaded files and the number of parsed files is patched to a small number where possible. However, most datasets require downloading one large file, which takes time. Hence removed from GitHub testing CI (by not naming it according to pytest convention).
'''

import unittest, tempfile, os, shutil, glob
from collections import defaultdict
from unittest import mock
from proteinshake.datasets import *

def download_mock(self):
    mock_data_path = os.path.dirname(os.path.realpath(__file__)) + '/mock_data'
    shutil.copytree(mock_data_path, f'{self.root}/raw/files', dirs_exist_ok=True)

def get_raw_files_mock(self):
    return glob.glob(f'{self.root}/raw/files/????.pdb')[:self.limit]

def get_id_from_filename_mock(self, filename):
    return filename[:4]

def scop_download_mock(self):
    download_mock(self)
    self.scop = self._parse_scop(f'{self.root}/raw/files/scop.txt')

def pli_download_mock(self):
    download_mock(self)
    self.index_data = self.parse_pdbbind_PL_index(f'{self.root}/raw/files/INDEX_refined_data.{self.version}')

def GODag_mock(*args, **kwargs):
    class DummyTerm:
        def __init__(self):
            self.namespace = 'molecular_function'
    factory = lambda: DummyTerm()
    return defaultdict(factory)

class TestDownload(unittest.TestCase):

    @mock.patch.object(ProteinLigandInterfaceDataset, 'download', pli_download_mock)
    @mock.patch.object(ProteinLigandInterfaceDataset, 'get_raw_files', get_raw_files_mock)
    def test_protein_ligand(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandInterfaceDataset(root=tmp, use_precomputed=False)
    
    @mock.patch.object(ProteinProteinInterfaceDataset, 'download', download_mock)
    @mock.patch.object(ProteinProteinInterfaceDataset, 'get_raw_files', get_raw_files_mock)
    def test_protein_protein(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinProteinInterfaceDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.RCSBDataset.limit', return_value=5)
    #@mock.patch.object(RCSBDataset, 'download', download_mock)
    #@mock.patch.object(RCSBDataset, 'get_raw_files', get_raw_files_mock)
    def test_rcsb(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            ds = RCSBDataset(root=tmp, use_precomputed=False)

    @mock.patch.object(GeneOntologyDataset, 'download', download_mock)
    @mock.patch.object(GeneOntologyDataset, 'get_raw_files', get_raw_files_mock)
    @mock.patch('proteinshake.datasets.gene_ontology.GODag', GODag_mock)
    def test_go(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = GeneOntologyDataset(root=tmp, use_precomputed=False)

    @mock.patch.object(EnzymeCommissionDataset, 'download', download_mock)
    @mock.patch.object(EnzymeCommissionDataset, 'get_raw_files', get_raw_files_mock)
    def test_ec(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = EnzymeCommissionDataset(root=tmp, use_precomputed=False)

    @mock.patch.object(ProteinFamilyDataset, 'download', download_mock)
    @mock.patch.object(ProteinFamilyDataset, 'get_raw_files', get_raw_files_mock)
    def test_pfam(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinFamilyDataset(root=tmp, use_precomputed=False)

    @mock.patch.object(TMAlignDataset, 'download', download_mock)
    @mock.patch.object(TMAlignDataset, 'get_raw_files', get_raw_files_mock)
    def test_tm(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = TMAlignDataset(root=tmp, use_precomputed=False)

    @mock.patch('proteinshake.datasets.AlphaFoldDataset.limit', return_value=5)
    #@mock.patch.object(AlphaFoldDataset, 'download', download_mock)
    #@mock.patch.object(AlphaFoldDataset, 'get_raw_files', get_raw_files_mock)
    #@mock.patch.object(AlphaFoldDataset, 'get_id_from_filename', get_id_from_filename_mock)
    def test_af(self, mock):
        with tempfile.TemporaryDirectory() as tmp:
            organism = 'methanocaldococcus jannaschii'
            ds = AlphaFoldDataset(root=tmp, organism=organism, use_precomputed=False)

    @mock.patch.object(SCOPDataset, 'download', scop_download_mock)
    @mock.patch.object(SCOPDataset, 'get_raw_files', get_raw_files_mock)
    def test_scop(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = SCOPDataset(root=tmp, use_precomputed=False)

    @mock.patch.object(ProteinLigandDecoysDataset, 'download', download_mock)
    @mock.patch.object(ProteinLigandDecoysDataset, 'get_raw_files', get_raw_files_mock)
    def test_decoys(self):
        with tempfile.TemporaryDirectory() as tmp:
            ds = ProteinLigandDecoysDataset(root=tmp, use_precomputed=False)

if __name__ == '__main__':
    unittest.main()
