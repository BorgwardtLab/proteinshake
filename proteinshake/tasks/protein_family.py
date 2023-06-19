from sklearn import metrics
import numpy as np

from proteinshake.datasets import ProteinFamilyDataset
from proteinshake.tasks import Task

class ProteinFamilyTask(Task):
    """ Predict the protein family classification of a protein structure. This is a protein-level multi-class prediction.
    """

    DatasetClass = ProteinFamilyDataset
    
    type = 'Multiclass Classification'
    input = 'Protein'
    output = 'Protein Family (Pfam)'
    default_metric = 'Precision'

    def compute_token_map(self):
        labels = {p['protein']['Pfam'][0] for p in self.dataset.proteins()}
        return {k: v for v,k in enumerate(sorted(labels))}

    def target(self, protein):
        return self.token_map[protein['protein']['Pfam'][0]]

    def evaluate(self, y_true, y_pred):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        y_true = np.array(y_true, dtype=int)
        y_pred = np.array(y_pred, dtype=int)
        return {
            'Precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'Recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'Accuracy': metrics.accuracy_score(y_true, y_pred),
        }
    
    def dummy(self):
        return np.random.choice(list(self.token_map.values()), len(self.test_targets))
