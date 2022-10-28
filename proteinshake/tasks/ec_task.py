from sklearn import metrics

from proteinshake.datasets import EnzymeCommissionDataset
from proteinshake.tasks import ShakeTask

class EnzymeCommissionTask(ShakeTask):
    """ Predict the enzyme commission classification of a protein structure. This is a protein-level multi-class prediction.

    """
    def __init__(self, root, ec_level=0, *args, **kwargs):
        dataset = EnzymeCommissionDataset(root=root)
        self.ec_level = ec_level

        super().__init__(dataset, *args, **kwargs)

    def compute_token_map(self):
        labels = {p['protein']['EC'].split(".")[self.ec_level] for p in self.proteins}
        self.token_map = {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def num_classes(self):
        return len(self.token_map)
    @property
    def task_type(self):
        return 'classification, multi-class'

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        return self.token_map[protein['protein']['EC'].split(".")[0]]

    def evaluate(self, pred, true):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        return {'precision': metrics.precision_score(pred, true, average='macro'),
                'recall': metrics.recall_score(pred, true, average='macro')}

if __name__ == "__main__":
    task = EnzymeCommissionTask(root='ec')
