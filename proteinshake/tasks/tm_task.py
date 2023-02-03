import itertools

import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split

from proteinshake.datasets import TMAlignDataset
from proteinshake.tasks import ShakeTask

class RetrieveTask(ShakeTask):
    """ Predict the structural similarity between two proteins. This is a pair-wise protein-level regression task.
    Ground truth is computed using the TMAlign software. Split indices are stored as tuples which contain two indices in
    the underlying dataset.
    """
    def __init__(self, root='data', *args, **kwargs):
        dataset = TMAlignDataset(root=root)
        super().__init__(dataset, *args, **kwargs)

    @property
    def task_type(self):
        return "regression"

    def compute_pairs(self, indices):
        com = list(itertools.combinations(range(len(indices)), 2))
        return np.array(indices)[com].tolist()

    def compute_random_split(self, *args, **kwargs):
        splits = super().compute_random_split(*args, **kwargs)
        return {k: self.compute_pairs(v) for k,v in splits.items()}

    def compute_cluster_split(self, *args, **kwargs):
        splits = super().compute_cluster_split(*args, **kwargs)
        return {k: self.compute_pairs(v) for k,v in splits.items()}

    def target(self, protein1, protein2=None):
        try:
            protein1, protein2 = protein1
        except:
            pass
        pdbid_1 = protein1['protein']['ID']
        pdbid_2 = protein2['protein']['ID']
        return self.dataset.tm_score[pdbid_1][pdbid_2]

    def evaluate(self, pred, true):
        return {'mse': metrics.mean_squared_error(pred, true)}
