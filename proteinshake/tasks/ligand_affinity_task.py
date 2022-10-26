from sklearn import metrics

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import ShakeTask

class LigandAffinityTask(ShakeTask):
    """ Predict the dissociation constant (kd) for a protein and a small molecule. This is a protein-level regression task.
    Small molecule ligand information is stored as ``dataset[i].smiles`` for a SMILES string, or as pre-computed molecular
    fingerprints ``dataset[i].fp_maccs``, ```dataset[i].fp_morgan_r2``.
    """
    def __init__(self, root, *args, **kwargs):
        dataset = ProteinLigandInterfaceDataset(root=root)
        super().__init__(dataset, *args, **kwargs)

    def target(self, idx):
        return self.proteins[idx]['kd']

    def evaluator(self, pred, true ):
        return {'mse': metrics.mean_squared_error(pred, true)}

if __name__ == "__main__":
    task = LigandAffinityTask(root='bob')
