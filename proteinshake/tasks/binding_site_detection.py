from sklearn import metrics

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import Task

class BindingSiteDetectionTask(Task):
    """ Identify the binding residues of a protein-small molecule binding site. This is a residue-level binary classification
    task.
    """
    
    DatasetClass = ProteinLigandInterfaceDataset

    @property
    def task_type(self):
        return "binary-classification"

    def target(self, protein):
        return protein['residue']['binding_site']

    def evaluate(self, pred, true):
        return {
                'accuracy': metrics.accuracy_score(true, pred),
                'mcc': metrics.matthews_corrcoef(true, pred),
                }
