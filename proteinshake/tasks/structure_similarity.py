from collections import defaultdict

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
    
    type = 'Regression'
    input = 'Protein and Protein'
    output = 'Local Distance Difference Test'
    default_metric = 'Spearman R'
    pairwise = True

    def target(self, protein_1, protein_2):
        pdbid_1 = protein_1['protein']['ID']
        pdbid_2 = protein_2['protein']['ID']
        return float(self.dataset.lddt(pdbid_1,pdbid_2))

    def evaluate(self, y_true, y_pred):
        return {
            'Mean Squared Error': metrics.mean_squared_error(y_true, y_pred),
            'Spearman R': spearmanr(y_true, y_pred)[0]
        }
    
    def dummy(self):
        return np.random.uniform(size=(len(self.test_targets),))
