import itertools

from scipy.stats import spearmanr
import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split

from proteinshake.datasets import TMAlignDataset
from proteinshake.tasks import Task

class StructureSimilarityTask(Task):
    """ Predict the structural similarity between two proteins. This is a pair-wise protein-level regression task.
    Ground truth is computed using the TMAlign software. Split indices are stored as tuples which contain two indices in
    the underlying dataset.
    """

    DatasetClass = TMAlignDataset

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def update_index(self):
        """ Transform to pairwise indexing """
        self.train_index = self.compute_pairs(self.train_index)
        self.val_index = self.compute_pairs(self.val_index)
        self.test_index = self.compute_pairs(self.test_index)
        
    def compute_targets(self):
        self.train_targets = np.array([self.target(*self.proteins[i]) for i in self.train_index])
        self.val_targets = np.array([self.target(*self.proteins[i]) for i in self.val_index])
        self.test_targets = np.array([self.target(*self.proteins[i]) for i in self.test_index])

    @property
    def task_type(self):
        return ('protein_pair', 'regression')

    def compute_pairs(self, index):
        combinations = np.array(list(itertools.combinations(range(len(index)), 2)), dtype=int)
        return index[combinations]

    def target(self, protein1, protein2):
        pdbid_1 = protein1['protein']['ID']
        pdbid_2 = protein2['protein']['ID']
        return self.dataset.lddt(pdbid_1,pdbid_2)

    def dummy_output(self):
        import random
        return [random.random() for _ in range(len(self.test_targets))]

    def evaluate(self, y_true, y_pred):
        return {
            'mse': metrics.mean_squared_error(y_true, y_pred),
            'r2':  spearmanr(y_true, y_pred)
        }
