import itertools
import numpy as np
from scipy.sparse import csr_matrix

from proteinshake.datasets import GeneOntologyDataset
from proteinshake.tasks import Task

class GeneOntologyTask(Task):
    """ Predict the Gene Ontology terms of a protein structure. This is a protein-level multi-label prediction.

    The prediction should be a n_samples x n_classes matrix, where the columns are ordered according to `self.classes`.
    If your model does not predict or handle a certain class, assign a zero value.
    """

    DatasetClass = GeneOntologyDataset
    
    type = 'Multilabel Classification'
    input = 'Protein'
    output = 'Gene Ontology Terms'
    default_metric = 'Fmax'

    def compute_token_map(self):
        labels = set(itertools.chain(*[p['protein']['molecular_function'] for p in self.dataset.proteins()]))
        return {k:v for v,k in enumerate(sorted(labels))}

    def target(self, protein):
        return [self.token_map[i] for i in protein['protein']['molecular_function']]
    
    def target_transform(self, index):
        transformed = np.zeros((len(index), len(self.token_map)), dtype=bool)
        targets = self.targets[index]
        for i,indices in enumerate(targets): transformed[i,indices] = True
        return transformed

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

    def fmax(self, y_true, y_pred):
        fmax = 0
        for t in np.linspace(0,1,21):
            prec, rec = self.precision(y_true, y_pred, t), self.recall(y_true, y_pred, t)
            if prec+rec == 0: continue
            f1 = (2 * prec * rec) / (prec + rec)
            fmax = max(fmax, f1)
        return fmax

    def evaluate(self, y_true, y_pred):
        y_pred = np.array(y_pred)
        return {
            'Fmax': self.fmax(y_true, y_pred),
        }
    
    def dummy(self):
        return np.random.rand(len(self.test_targets), len(self.token_map))
