from .task import Task
from .enzyme_class import EnzymeClassTask
from .ligand_affinity import LigandAffinityTask
from .binding_site__detection import BindingSiteDetectionTask
from .structure_similarity import StructureSimilarityTask
from .structural_class import StructuralClassTask

classes = ['Task',
           'EnzymeClassTask',
           'LigandAffinityTask',
           'StructureSimilarityTask',
           'BindingSiteDetectionTask',
           'StructuralClassTask'
           ]

__all__ = classes
