from sklearn import metrics
import numpy as np

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import Task

class LigandAffinityTask(Task):
    """ Predict the dissociation constant (Kd) for a protein and a small molecule. This is a protein-level regression task.
    Small molecule ligand information is stored as ``dataset[i].smiles`` for a SMILES string, or as pre-computed molecular
    fingerprints ``dataset[i].fp_maccs``, ```dataset[i].fp_morgan_r2``.
    """
    
    DatasetClass = ProteinLigandInterfaceDataset
    
    type = 'Regression'
    input = 'Protein and Molecule'
    output = 'Dissociation Constant Kd'
    default_metric = 'R2'

    def target(self, protein_dict):
        return protein_dict['protein']['neglog_aff']

    def evaluate(self, y_true, y_pred):
        return {
            'Mean Squared Error': metrics.mean_squared_error(y_true, y_pred),
            'R2': metrics.r2_score(y_true, y_pred)
        }
    
    def dummy(self):
        return np.random.uniform(len(self.test_targets))
