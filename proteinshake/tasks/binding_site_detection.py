from sklearn import metrics

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import Task

class BindingSiteDetectionTask(Task):
    """ Identify the binding residues of a protein-small molecule binding site. This is a residue-level binary classification task.

    """
    
    DatasetClass = ProteinLigandInterfaceDataset

    @property
    def task_type(self):
        return ('residue', 'binary')

    def dummy_output(self):
        import random
        return [random.randint(0, 1) for p in self.test_targets]

    def target(self, protein):
        return protein['residue']['binding_site']

    def compute_targets(self):
        # compute targets (e.g. for scaling)
        self.train_targets = [p for i in self.train_index for p in self.target(self.proteins[i])]
        self.val_targets = [p for i in self.val_index for p in self.target(self.proteins[i])]
        self.test_targets = [p for i in self.test_index for p in self.target(self.proteins[i])]

    def evaluate(self, y_true, y_pred):
        return {
            'accuracy': metrics.accuracy_score(y_true, y_pred),
            'mcc': metrics.matthews_corrcoef(y_true, y_pred),
        }
