from sklearn import metrics
from functools import cached_property
import numpy as np

from proteinshake.datasets import EnzymeCommissionDataset
from proteinshake.tasks import Task

class EnzymeClassTask(Task):
    """ Predict the enzyme commission classification of a protein structure. This is a protein-level multi-class prediction.

    """

    DatasetClass = EnzymeCommissionDataset
    
    type = 'Multiclass Classification'
    input = 'Protein'
    output = 'Enzyme Commission Level 1'
    
    def __init__(self, ec_level=0, *args, **kwargs):
        self.ec_level = ec_level
        super().__init__(*args, **kwargs)

    @cached_property
    def token_map(self):
        labels = {p['protein']['EC'].split(".")[self.ec_level] for p in self.proteins}
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    def dummy_output(self):
        import random
        tokens = list(self.token_map.values())
        return [random.choice(tokens) for _ in range(len(self.test_targets))]

    @property
    def num_classes(self):
        return len(self.token_map)

    @property
    def task_type(self):
        return ('protein', 'multi_class')

    @property
    def task_in(self):
        return ('protein')

    @property
    def task_out(self):
        return ('multi_class')

    @property
    def target_dim(self):
        return (len(self.token_map.values()))

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        return self.token_map[protein['protein']['EC'].split(".")[0]]

    @property
    def default_metric(self):
        return 'precision'

    def evaluate(self, y_true, y_pred):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        y_true = np.array(y_true, dtype=int)
        y_pred = np.array(y_pred, dtype=int)
        return {
            'precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'accuracy': metrics.accuracy_score(y_true, y_pred),
            #'AUROC': metrics.roc_auc_score(y_true, y_pred, average='macro', multi_class='ovo'),
        }
