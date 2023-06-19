from sklearn import metrics
import numpy as np

from proteinshake.datasets import EnzymeCommissionDataset
from proteinshake.tasks import Task

class EnzymeClassTask(Task):
    """ Predict the enzyme commission classification of a protein structure.
    """

    DatasetClass = EnzymeCommissionDataset
    
    type = 'Multiclass Classification'
    input = 'Protein'
    output = 'Enzyme Commission Level 1'
    default_metric = 'Precision'

    def compute_token_map(self):
        unique_labels = {p['protein']['EC'].split(".")[0] for p in self.dataset.proteins()}
        return {k:v for v,k in enumerate(sorted(unique_labels))}

    def target(self, protein_dict):
        return self.token_map[protein_dict['protein']['EC'].split(".")[0]]

    def evaluate(self, y_true, y_pred):
        y_true = np.array(y_true, dtype=int)
        y_pred = np.array(y_pred, dtype=int)
        return {
            'Precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'Recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'Accuracy': metrics.accuracy_score(y_true, y_pred),
        }
    
    def dummy(self):
        return np.random.choice(list(self.token_map.values()), len(self.test_targets))
