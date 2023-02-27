from .task import Task
from .enzyme_class import EnzymeClassTask
from .ligand_affinity import LigandAffinityTask
from .binding_site_detection import BindingSiteDetectionTask
from .structure_similarity import StructureSimilarityTask
from .structural_class import StructuralClassTask
from .gene_ontology import GeneOntologyTask

classes = ['Task',
           'EnzymeClassTask',
           'LigandAffinityTask',
           'StructureSimilarityTask',
           'BindingSiteDetectionTask',
           'StructuralClassTask',
           'GeneOntologyTask',
           ]

__all__ = classes
