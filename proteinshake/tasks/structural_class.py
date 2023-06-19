from sklearn import metrics
import numpy as np

from proteinshake.datasets import SCOPDataset
from proteinshake.tasks import Task

class StructuralClassTask(Task):
    """ Predict the SCOP class of a protein structure. This is a protein-level multi-class prediction.
    """

    DatasetClass = SCOPDataset
    
    type = 'Multiclass Classification'
    input = 'Protein'
    output = 'SCOP Class'
    default_metric = 'Accuracy'
    
    def __init__(self, scop_level='SCOP-FA', *args, **kwargs):
        self.scop_level = scop_level
        super().__init__(*args, **kwargs)

    def compute_token_map(self):
        labels = {p['protein'][self.scop_level] for p in self.dataset.proteins()}
        return {k:v for v,k in enumerate(sorted(list(labels)))}

    def target(self, protein):
        return self.token_map[protein['protein'][self.scop_level]]

    def evaluate(self, y_true, y_pred):
        y_true = np.array(y_true, dtype=int)
        y_pred = np.array(y_pred, dtype=int)
        return {
            'Precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'Recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'Accuracy': metrics.accuracy_score(y_true, y_pred),
        }
    
    def dummy(self):
        return np.random.choice(self.token_map.values(), size=len(self.test_targets))
