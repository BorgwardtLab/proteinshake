from sklearn import metrics

from proteinshake.datasets import EnzymeCommissionDataset
from proteinshake.tasks import ShakeTask

class EnzymeCommissionTask(ShakeTask):
    def __init__(self, root, *args, **kwargs):
        dataset = EnzymeCommissionDataset(root)

        super().__init__(dataset, *args, **kwargs)

    def compute_token_map(self):
        proteins, _ = self.dataset.proteins()
        labels = {self.label_process(p['protein'][self.target]) for p in proteins}
        self.token_map = {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def target(self, idx):
        return self.token_map[self.proteins[idx]['EC'].split(".")[0]]

    def evaluate(self, pred, true):
        return {'precision': metrics.precision_score(pred, true)}

    @property
    def level(self):
        return 'protein'
