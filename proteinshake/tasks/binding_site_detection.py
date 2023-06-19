from sklearn import metrics
import numpy as np

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import Task

class BindingSiteDetectionTask(Task):
    """ Identify the binding residues of a protein-small molecule binding site. This is a residue-level binary classification task.

    """
    
    DatasetClass = ProteinLigandInterfaceDataset
    
    type = 'Binary Classification'
    input = 'Protein'
    output = 'Small Molecule Binding Residues'
    default_metric = 'MCC'
    level = 'Residue'

    def target(self, protein_dict):
        binding_sites = protein_dict['residue']['binding_site']
        return np.arange(len(binding_sites))[binding_sites].tolist()

    def evaluate(self, y_true, y_pred):
        return {
            'Accuracy': metrics.accuracy_score(y_true.flatten(), y_pred.flatten()),
            'MCC': metrics.matthews_corrcoef(y_true.flatten(), y_pred.flatten()),
        }
    
    def dummy(self):
        return [np.random.choice([0,1], size=L) for L in self.sizes[self.test_index]]
