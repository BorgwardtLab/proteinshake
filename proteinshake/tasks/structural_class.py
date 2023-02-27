from sklearn import metrics

from proteinshake.datasets import SCOPDataset
from proteinshake.tasks import Task

class StructuralClassTask(Task):
    """ Predict the SCOP class of a protein structure. This is a protein-level multi-class prediction.

    """

    DatasetClass = SCOPDataset
    
    def __init__(self, scop_level='SCOP-FA', *args, **kwargs):
        self.scop_level = scop_level
        super().__init__(*args, **kwargs)

    def compute_token_map(self):
        labels = {p['protein'][self.scop_level] for p in self.proteins()}
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def num_classes(self):
        return len(self.token_map)

    @property
    def task_type(self):
        return 'classification, multi-class'

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        return self.token_map[protein['protein'][self.scop_level]]

    def evaluate(self, y_pred):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        return {
            'precision': metrics.precision_score(self.test_targets, y_pred, average='macro', zero_division=0),
            'recall': metrics.recall_score(self.test_targets, y_pred, average='macro', zero_division=0),
            'accuracy': metrics.accuracy_score(self.test_targets, y_pred),
            #'AUROC': metrics.roc_auc_score(self.test_targets, y_pred, average='macro', multi_class='ovo'),
        }
