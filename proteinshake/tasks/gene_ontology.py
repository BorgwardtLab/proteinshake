import itertools
import numpy as np
from sklearn import metrics
from functools import cached_property

from proteinshake.datasets import GeneOntologyDataset
from proteinshake.tasks import Task

class GeneOntologyTask(Task):
    """ Predict the Gene Ontology terms describing the functional roles of a given protein in the cell. This is a protein-level multi-label prediction.

    The prediction should be a n_samples x n_classes matrix, where the columns are ordered according to `self.classes`.
    If your model does not predict or handle a certain class, assign a zero value.

    .. admonition:: Task Summary

        * **Input:** one protein
        * **Output:** n_classes gene ontology terms 
        * **Evaluation:** Fmax (Radivojac, Predrag, et al. "A large-scale evaluation of computational protein function prediction." Nature methods 10.3 (2013): 221-227.)


    """

    DatasetClass = GeneOntologyDataset
    
    type = 'Multilabel Classification'
    input = 'Protein'
    output = 'Gene Ontology Terms'
    
    def __init__(self, branch='molecular_function', *args, **kwargs):
        self.branch = branch
        super().__init__(*args, **kwargs)

    @cached_property
    def token_map(self):
        labels = set(itertools.chain(*[p['protein'][self.branch] for p in self.proteins]))
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def num_classes(self):
        return len(self.token_map)
    
    @property
    def classes(self):
        return list(self.token_map.keys())

    @property
    def task_in(self):
        return ('protein')

    @property
    def task_type(self):
        return ('protein', 'multi_label')

    @property
    def task_out(self):
        return ('multi_label')

    @property
    def target_dim(self):
        return (len(self.token_map.values()))

    @property
    def num_features(self):
        return 20

    def target(self, protein):
        tokens = [self.token_map[i] for i in protein['protein'][self.branch]]
        target = np.zeros_like(self.classes, dtype=bool)
        target[tokens] = True
        return target

    def precision(self, y_true, y_pred, threshold):
        mt = (y_pred.max(axis=1) >= threshold).sum()
        if mt == 0: return 0.0
        y_pred = y_pred >= threshold
        nom = np.logical_and(y_true, y_pred).sum(axis=1).astype(np.float32)
        denom = y_pred.sum(axis=1).astype(np.float32)
        return 1/mt * np.divide(nom, denom, out=np.zeros_like(nom), where=denom!=0).sum()

    def recall(self, y_true, y_pred, threshold):
        ne = y_true.shape[0]
        if ne == 0: return 0.0
        y_pred = y_pred >= threshold
        nom = np.logical_and(y_true, y_pred).sum(axis=1).astype(np.float32)
        denom = y_true.sum(axis=1).astype(np.float32)
        return 1/ne * np.divide(nom, denom, out=np.zeros_like(nom), where=denom!=0).sum()
    
    def remaining_uncertainty(self, y_true, y_pred, threshold):
        pass

    def missing_information(self, y_true, y_pred, threshold):
        pass

    def fmax(self, y_true, y_pred):
        fmax = 0
        for t in np.linspace(0,1,21):
            prec, rec = self.precision(y_true, y_pred, t), self.recall(y_true, y_pred, t)
            if prec+rec == 0: continue
            f1 = (2 * prec * rec) / (prec + rec)
            fmax = max(fmax, f1)
        return fmax
    
    def smin(self, y_pred):
        return min([
            np.sqrt(
                self.remaining_uncertainty(y_pred, t) ** 2
                + self.missing_information(y_pred, t) ** 2
            )
            for t in np.linspace(0,1,21)
        ])

    def dummy_output(self):
        return np.random.rand(len(self.test_index), len(self.token_map.keys()))

    @property
    def default_metric(self):
        return 'Fmax'

    def evaluate(self, y_true, y_pred):
        y_true, y_pred = np.array(y_true), np.array(y_pred)
        return {
            'Fmax': self.fmax(y_true, y_pred),
            #'Smin': self.smin(y_true, y_pred),
        }
