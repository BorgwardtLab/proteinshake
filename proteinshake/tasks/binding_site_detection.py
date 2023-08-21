from sklearn import metrics

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import Task

class BindingSiteDetectionTask(Task):
    """ Identify the binding residues (binding pocket) of a protein-small molecule binding site. An important step in drug discovery for proteins is to find 
    potential cavities where small molecules can bind the protein based on the whole protein's structure. 
 Pocket atoms/residues
    taken directly from PDBBind annotations.

    .. admonition:: Task Summary 

        * **Input:** one protein
        * **Output:** binary label for each atom/residue
        * **Evaluation:** Matthews Correlation Coefficient (*Gallo Cassarino, Tiziano, Lorenza Bordoli, and Torsten Schwede. "Assessment of ligand binding site predictions in CASP10." Proteins: Structure, Function, and Bioinformatics 82 (2014): 154-163.*)


    """
    
    DatasetClass = ProteinLigandInterfaceDataset
    
    type = 'Binary Classification'
    input = 'Residue'
    output = 'Small Molecule Binding Residues'
    
    @property
    def num_classes(self):
        return 2

    @property
    def task_in(self):
        return ('residue')

    @property
    def task_type(self):
        return ('residue', 'binary')

    @property
    def task_out(self):
        return ('binary')
    
    @property
    def target_dim(self):
        return (1)

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

    @property
    def default_metric(self):
        return 'mcc'

    def evaluate(self, y_true, y_pred):
        return {
            'accuracy': metrics.accuracy_score(y_true, y_pred),
            'mcc': metrics.matthews_corrcoef(y_true, y_pred),
        }
