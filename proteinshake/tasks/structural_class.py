from sklearn import metrics
from functools import cached_property
import numpy as np

from proteinshake.datasets import SCOPDataset
from proteinshake.tasks import Task

class StructuralClassTask(Task):
    """ Predict the SCOP class of a protein structure. SCOP labels proteins according to a hierarchy of structural and evolutionary information. The top level of the hierarchy ``SCOP_FA``, you can customize the task to use a different level setting ``scop_level`` to ``SCOP_{level}``, where level is any of TP=protein type, CL=protein class, CF=fold, SF=superfamily, FA=family. This is a protein-level multi-class prediction.

    .. admonition:: Task Summary 

        * **Input:** one protein
        * **Output:** SCOP class (3042 classes)
        * **Evaluation:** Accuracy (custom task)

    """

    DatasetClass = SCOPDataset
    
    type = 'Multiclass Classification'
    input = 'Protein'
    output = 'SCOP Class'
    
    def __init__(self, scop_level='SCOP-FA', *args, **kwargs):
        self.scop_level = scop_level
        super().__init__(*args, **kwargs)
        
    @property
    def num_classes(self):
        return len(self.token_map)

    @cached_property
    def token_map(self):
        labels = {p['protein'][self.scop_level] for p in self.proteins}
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def target_dim(self):
        return len(self.token_map)

    def dummy_output(self):
        import random
        tokens = list(self.token_map.values())
        return [random.choice(tokens) for _ in range(len(self.test_targets))]

    @property
    def task_in(self):
        return ('protein')

    @property
    def task_type(self):
        return ('protein', 'multi-class')

    @property
    def task_out(self):
        return ('multi_class')

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        return self.token_map[protein['protein'][self.scop_level]]

    @property
    def default_metric(self):
        return 'accuracy'

    def evaluate(self, y_true, y_pred):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        y_true = np.array(y_true, dtype=int)
        y_pred = np.array(y_pred, dtype=int)
        return {
            'precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'accuracy': metrics.accuracy_score(y_true, y_pred),
            #'AUROC': metrics.roc_auc_score(self.test_targets, y_pred, average='macro', multi_class='ovo'),
        }
