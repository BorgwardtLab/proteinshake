from .task import ShakeTask
from .ec_task import EnzymeCommissionTask
from .ligand_affinity_task import LigandAffinityTask
from .binding_site_task import BindingSitePredictionTask
from .tm_task import RetrieveTask

classes = ['ShakeTask',
           'EnzymeCommissionTask',
           'LigandAffinityTask',
           'RetrieveTask',
           'BindingSitePredictionTask'
           ]

__all__ = classes
