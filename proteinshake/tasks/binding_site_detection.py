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
        return np.arange(len(binding_sites))[binding_sites].tolist(), len(binding_sites)
    
    def target_transform(self, target):
        residue_indices, size = target
        transformed_target = np.zeros(size, dtype=bool)
        transformed_target[residue_indices] = True
        return transformed_target

    def evaluate(self, y_true, y_pred):
        print(np.hstack(y_true))
        print(np.hstack(y_pred))
        return {
            'Accuracy': metrics.accuracy_score(np.hstack(y_true), np.hstack(y_pred)),
            'MCC': metrics.matthews_corrcoef(np.hstack(y_true), np.hstack(y_pred)),
        }
    
    @property
    def y_dummy(self):
        return np.array([np.random.choice([0,1], size=len(x)).astype(bool) for x in self.y_test], dtype=object)
