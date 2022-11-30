from .dataset import Dataset
from .rcsb import RCSBDataset
from .enzyme_commission import EnzymeCommissionDataset
from .gene_ontology import GeneOntologyDataset
from .protein_protein_interface import ProteinProteinInterfaceDataset
from .protein_ligand_interface import ProteinLigandInterfaceDataset
from .pfam import PfamDataset
from .tm_align import TMAlignDataset
from .alphafold import AlphaFoldDataset
from .scop import SCOPDataset

__all__ = [
    'Dataset',
    'RCSBDataset',
    'EnzymeCommissionDataset',
    'GeneOntologyDataset',
    'ProteinProteinInterfaceDataset',
    'ProteinLigandInterfaceDataset',
    'PfamDataset',
    'TMAlignDataset',
    'AlphaFoldDataset',
    'SCOPDataset'
    ]

classes = __all__
