from sklearn import metrics

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import Task

class LigandAffinityTask(Task):
    """ Predict the dissociation constant (Kd) for a protein and a small molecule. Accurate estimates of the binding strength between a protein and a small molecule ligand is a crucial step in understanding protein function regulation and in efficiently searching the massive space of small molecules new therapies. 
    Small molecule ligand information is stored as ``dataset[i].smiles`` for a SMILES string, or as pre-computed molecular
    fingerprints ``dataset[i].fp_maccs``, ```dataset[i].fp_morgan_r2``.

    .. admonition:: Task Summary 

        * **Input:** One protein and one ligand SMILES string
        * **Output:** predicted dissociation constant (scalar)
        * **Evaluation:** R2 score (*Stepniewska-Dziubinska, Marta M., Piotr Zielenkiewicz, and Pawel Siedlecki. "Development and evaluation of a deep learning model for protein–ligand binding affinity prediction." Bioinformatics 34.21 (2018): 3666-3674.*) 

    """
    
    DatasetClass = ProteinLigandInterfaceDataset
    
    type = 'Regression'
    input = 'Protein and Molecule'
    output = 'Dissociation Constant Kd'

    @property
    def task_in(self):
        return ('protein', 'molecule')

    @property
    def task_type(self):
        return ('protein', 'regression')

    @property
    def task_out(self):
        return ('regression')

    @property
    def out_dim(self):
        return (1)

    def dummy_output(self):
        import random
        return [random.random() for _ in range(len(self.test_index))] 

    def target(self, protein):
        return protein['protein']['neglog_aff']

    @property
    def default_metric(self):
        return 'r2'

    def evaluate(self, y_true, y_pred):
        return {
            'mse': metrics.mean_squared_error(y_true, y_pred),
            'r2': metrics.r2_score(y_true, y_pred)
        }
