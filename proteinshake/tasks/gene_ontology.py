import itertools
import numpy as np
from sklearn import metrics

from proteinshake.datasets import GeneOntologyDataset
from proteinshake.tasks import Task

class GeneOntologyTask(Task):
    """ Predict the Gene Ontology terms of a protein structure. This is a protein-level multi-label prediction.

    The prediction should be a n_samples x n_classes matrix, where the columns are ordered according to `self.classes`.
    If your model does not predict or handle a certain class, assign a zero value.
    """

    DatasetClass = GeneOntologyDataset
    
    def __init__(self, branch='molecular_function', *args, **kwargs):
        self.branch = branch
        super().__init__(*args, **kwargs)

    def compute_token_map(self):
        labels = set(itertools.chain(*[p['protein'][self.branch] for p in self.proteins]))
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def num_classes(self):
        return len(self.token_map)
    
    @property
    def classes(self):
        return list(self.token_map.keys())

    @property
    def task_type(self):
        return 'classification, multi-label'

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        return [term in protein['protein'][self.branch] for term in self.classes]

    def precision(self, y_pred, threshold):
        y_true = np.copy(self.test_targets)[y_pred.max(axis=1) >= threshold]
        y_pred = y_pred[y_pred.max(axis=1) >= threshold] >= threshold
        mt = y_pred.shape[0]
        return 1/mt * (
            np.logical_and(y_true, y_pred).sum(axis=1)
            / y_pred.sum(axis=1)
            ).sum()


    def recall(self, y_pred, threshold):
        y_true = np.copy(self.test_targets)
        y_pred = y_pred >= threshold
        ne = y_true.shape[0]
        return 1/ne * (
            np.logical_and(y_true, y_pred).sum(axis=1)
            / y_true.sum(axis=1)
            ).sum()
    
    def remaining_uncertainty(self, y_pred, threshold):
        pass

    def missing_information(self, y_pred, threshold):
        pass

    def fmax(self, y_pred):
        return max([
            2 * self.precision(y_pred, t) * self.recall(y_pred, t)
            / (self.precision(y_pred, t) + self.recall(y_pred, t))
            for t in np.linspace(0,1,21)
        ])
    
    def smin(self, y_pred):
        return min([
            np.sqrt(
                self.remaining_uncertainty(y_pred, t) ** 2
                + self.missing_information(y_pred, t) ** 2
            )
            for t in np.linspace(0,1,21)
        ])

    def evaluate(self, y_pred):
        return {
            'Fmax': self.fmax(y_pred),
            #'Smin': self.smin(y_pred),
        }
