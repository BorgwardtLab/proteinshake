from sklearn import metrics
from functools import cached_property

from proteinshake.datasets import ProteinFamilyDataset
from proteinshake.tasks import Task

class ProteinFamilyTask(Task):
    """ Predict the protein family classification of a protein structure. This is a protein-level multi-class prediction.

    """

    DatasetClass = ProteinFamilyDataset
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    @cached_property
    def token_map(self):
        # Pfam': ['Fis1 N-terminal tetratricopeptide repeat (Fis1_TPR_N)', 'Fis1 C-terminal tetratricopeptide repeat (Fis1_TPR_C)'], 
        labels = {p['protein']['Pfam'][0] for p in self.proteins}
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def num_classes(self):
        return len(self.token_map)

    def dummy_output(self):
        import random
        tokens = list(self.token_map.values())
        return [random.choice(tokens) for _ in range(len(self.test_targets))]

    @property
    def task_type(self):
        return ('protein', 'multi_class')

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        return self.token_map[protein['protein']['Pfam'][0]]

    def evaluate(self, y_true, y_pred):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        return {
            'precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'accuracy': metrics.accuracy_score(y_true, y_pred),
            #'AUROC': metrics.roc_auc_score(y_true, y_pred, average='macro', multi_class='ovo'),
        }
