from proteinshake.tasks.task import ShakeTask
from sklearn import metrics

class LigandAffinityTask(ShakeTask):
    def __init__(self, dataset, *args, **kwargs):
        super().__init__(dataset, *args, **kwargs)

    @property
    def target(self, idx):
        return self.proteins[idx]['kd']

    def evaluator(self):
        return metrics.mean_squared_error
