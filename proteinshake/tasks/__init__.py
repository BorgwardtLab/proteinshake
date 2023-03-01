from .task import Task
from .enzyme_class import EnzymeClassTask
from .ligand_affinity import LigandAffinityTask
from .binding_site_detection import BindingSiteDetectionTask
from .structure_similarity import StructureSimilarityTask
from .ppi_task import ProteinProteinInterfaceTask
from .tm_search import StructureSearchTask
from .structural_class import StructuralClassTask

classes = ['Task',
           'EnzymeClassTask',
           'LigandAffinityTask',
           'BindingSiteDetectionTask',
           'ProteinProteinInterfaceTask',
           'StructureSimilarityTask',
           'StructureSearchTask',
           'StructuralClassTask'
           ]

__all__ = classes
