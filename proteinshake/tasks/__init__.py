from .task import ShakeTask
from .ec_task import EnzymeCommissionTask
from .ligand_affinity_task import LigandAffinityTask
from .binding_site_task import BindingSitePredictionTask
from .ppi_task import ProteinProteinInterfaceTask
from .tm_task import RetrieveTask
from .tm_search import StructureSearchTask
from .scop_task import SCOPTask

classes = ['ShakeTask',
           'EnzymeCommissionTask',
           'LigandAffinityTask',
           'ProteinProteinInterfaceTask',
           'RetrieveTask',
           'BindingSitePredictionTask',
           'StructureSearchTask',
           'SCOPTask'
           ]

__all__ = classes
