from .task import ShakeTask
from .ec_task import EnzymeCommissionTask
from .ligand_affinity_task import LigandAffinityTask
from .tm_task import RetrieveTask

classes = ['ShakeTask',
           'EnzymeCommissionTask',
           'LigandAffinityTask',
           'RetrieveTask'
           ]

__all__ = classes
