from sklearn import metrics

from proteinshake.datasets import EnzymeCommissionDataset
from proteinshake.tasks import ShakeTask

class EnzymeCommissionTask(ShakeTask):
    """ Predict the enzyme commission classification of a protein structure. This is a protein-level multi-class prediction.

    """
    def __init__(self, root, *args, **kwargs):
        dataset = EnzymeCommissionDataset(root=root)

        super().__init__(dataset, *args, **kwargs)

    def compute_token_map(self):
        labels = {p['protein']['EC'].split(".")[0] for p in self.proteins}
        self.token_map = {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def num_classes(self):
        return len(self.token_map)
    @property
    def task_type(self):
        return 'classification, multi-class'

    def target(self, idx):
        return self.token_map[self.proteins[idx]['protein']['EC'].split(".")[0]]

    def evaluate(self, pred, true):
        return {'precision': metrics.precision_score(pred, true, average='macro')}
