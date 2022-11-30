from sklearn import metrics

from proteinshake.datasets import SCOPDataset
from proteinshake.tasks import ShakeTask

class SCOPTask(ShakeTask):
    """ Predict the enzyme commission classification of a protein structure. This is a protein-level multi-class prediction.

    """
    def __init__(self, root, scop_level='SCOP-FA', *args, **kwargs):
        dataset = SCOPDataset(root=root)
        self.scop_level = scop_level

        super().__init__(dataset, *args, **kwargs)

    def compute_token_map(self):
        labels = {p['protein'][self.scop_level] for p in self.proteins}
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
        return self.token_map[protein['protein'][self.scop_level]]

    def evaluate(self, y_true, y_pred):
        """ Using metrics from https://doi.org/10.1073/pnas.1821905116 """
        return {
            'precision': metrics.precision_score(y_true, y_pred, average='macro', zero_division=0),
            'recall': metrics.recall_score(y_true, y_pred, average='macro', zero_division=0),
            'accuracy': metrics.accuracy_score(y_true, y_pred),
            #'AUROC': metrics.roc_auc_score(y_true, y_pred, average='macro', multi_class='ovo'),
        }

if __name__ == "__main__":
    task = SCOPTask(root='scop')
