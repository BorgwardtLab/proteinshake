import itertools
from collections import defaultdict

from scipy.stats import spearmanr
import numpy as np
from sklearn import metrics
from sklearn.model_selection import train_test_split

from proteinshake.datasets import TMAlignDataset
from proteinshake.tasks import ShakeTask

class StructureSearchTask(ShakeTask):
    """ Retrieve similar proteins to a query.
    Evaluation is cast in the setting of recommender systems where we wish
    to retrieve 'relevant' documents from a large pool of documents.
    Here, a protein is a document and the relevant ones are all proteins
    with a minimum similarity to the query protein.

    """
    def __init__(self, min_sim=0.8, root='data', *args, **kwargs):
        dataset = TMAlignDataset(root=root)
        self.dataset = dataset
        self.min_sim = min_sim
        self.targets = self.compute_targets()

        super().__init__(dataset, *args, **kwargs)

    @property
    def task_type(self):
        return "retrieval"

    def target(self, protein):
        """ The target for a protein is a list of  proteins deemed 'relevant'
        according to `self.min_sim`.
        """
        return self.targets[protein['protein']['ID']]

    def compute_targets(self):
        """ Precompute the set of similar proteins for each query """
        targets = {}
        for q, candidates in self.dataset.tm_score.items():
            targets[q] = [c for c, sim in candidates.items() if sim >= self.min_sim]
        return targets

    def _precision_at_k(self, pred, targets, k):
        return len(set(pred[:k]).intersection(set(targets))) / len(pred)

    def _recall_at_k(self, pred, targets, k):
        return len(set(pred[:k]).intersection(set(targets))) / len(targets)

    def evaluate(self, pred, k=5):
        """ Retrieval metrics.

        Arguments
        -----------
        pred:
            List of indices of items in the dataset for which predictions were made.
        """
        results = defaultdict(list)
        for query, preds in zip(self.proteins[self.test_index], pred):
            targets = self.target(query)
            results['precision_at_k'].append(self._precision_at_k(preds, targets, k))
            results['recall_at_k'].append(self._recall_at_k(preds, targets, k))

        return {k: np.mean(v) for k, v in results.items()}

if __name__ == "__main__":
    ta = StructureSearchTask(use_precomputed=False).to_graph(eps=8).pyg()
    k = 5
    # generate random predictions
    ids = [p['protein']['ID'] for p in ta.proteins]
    pred = [sample(ids, k) for _ in range(len(ids))]
    metrics = ta.evaluate(pred, k=k)
    print(metrics)
