from sklearn import metrics

from proteinshake.datasets import ProteinLigandInterfaceDataset
from proteinshake.tasks import ShakeTask

class LigandAffinityTask(ShakeTask):
    def __init__(self, root, *args, **kwargs):
        dataset = ProteinLigandInterfaceDataset(root=root)
        super().__init__(dataset, *args, **kwargs)

    @property
    def target(self, idx):
        return self.proteins[idx]['kd']

    def evaluator(self, pred, true ):
        return {'mse': metrics.mean_squared_error(pred, true)}

if __name__ == "__main__":
    task = LigandAffinityTask(root='bob')
