import os
import json
import itertools

import numpy as np
from sklearn.model_selection import train_test_split

from proteinshake.utils import tmalign_wrapper, cdhit_wrapper, download_url, save, load

class Task:
    """ Base class for task-related utilities.
    This class wraps a proteinshake dataset and exposes split indices,
    integer-coded labels for classification tasks, and an evaluator function.
    Users should use this class to build their own dataloaders for training
    and evaluating models.

    Sample usage (assuming you have a model and dataloader in the namespace):

     .. code-block:: python

        >>> from proteinshake.tasks import EnzymeCommissionTask
        >>> task = EnzymeCommissionTask()
        >>> y_pred = model(task.train)
        >>> task.evaluate(y_pred)
        ... {'roc_auc_score': 0.7}


    Parameters
    ----------
    dataset: pytorch.datasets.Dataset
        Dataset to use for this task.
    split: str, default='random'
        How to split the data. Can be 'random', 'sequence', or 'structure'.
    random_state: int, default=42
        Random seed for reproducible splitting.
    train_ratio: float, default=0.80
        Fraction of dataset to use for training.
    val_ratio: float, default=0.10
        Fraction of dataset to use for validation.
    test_ratio: float, default=0.10
        Fraction of dataset to use for testing.
    cace_dir: str, default='.proteinshake'
        Directory where we store the result of computing splits and tokenizing.
    use_precomputed: bool, default=True
    """

    DatasetClass = None

    def __init__(self,
                 root               = 'data',
                 split              = 'random',
                 random_state       = 42,
                 train_ratio        = 0.80,
                 val_ratio          = 0.10,
                 test_ratio         = 0.10,
                 use_precomputed    = True,
                 **kwargs
                ):
        self.dataset = self.DatasetClass(root=root)
        self.proteins, self.size = self.dataset.proteins()
        class Proteins():
            def __init__(self, proteins):
                self.proteins = list(proteins)
            def __getitem__(self, idx):
                try:
                    idx = int(idx)
                except:
                    return [self.__getitem__(i) for i in idx]
                if idx >= len(self.proteins):
                    raise StopIteration
                return self.proteins[idx]
        self.proteins = Proteins(self.proteins)
        self.name = self.__class__.__name__

        # load task info
        path = f'{self.dataset.root}/{self.name}.json.gz'
        if not os.path.exists(path):
            if use_precomputed:
                self.download_precomputed()
            else:
                self.compute_splits(train_ratio, val_ratio, test_ratio, random_state)
        info = load(path)

        # load split indices
        self.train_index = np.array(info['splits'][split]['train'])
        self.val_index = np.array(info['splits'][split]['val'])
        self.test_index = np.array(info['splits'][split]['test'])
        self.token_map = info['token_map']

        # compute targets (e.g. for scaling)
        self.train_targets = [self.target(self.proteins[i]) for i in self.train_index]
        self.val_targets = [self.target(self.proteins[i]) for i in self.val_index]
        self.test_targets = [self.target(self.proteins[i]) for i in self.test_index]

    def download_precomputed(self):
        download_url(f'{self.dataset.repository_url}/{self.name}.json.gz', f'{self.dataset.root}')

    def compute_splits(self, *args, **kwargs):
        print('Computing splits...')
        info = {
            'splits': {
                'random': self.compute_random_split(*args, **kwargs),
                'sequence': self.compute_cluster_split('sequence', *args, **kwargs),
                #'structure': self.compute_cluster_split('structure', *args, **kwargs)
            },
            'token_map': self.compute_token_map()
        }
        save(info, f'{self.dataset.root}/{self.name}.json.gz')

    def compute_random_split(self, train_ratio, val_ratio, test_ratio, random_state):
        inds = range(self.size)
        train, test = train_test_split(inds, test_size=1-train_ratio, random_state=random_state)
        val, test = train_test_split(test, test_size=test_ratio/(test_ratio+val_ratio), random_state=random_state)
        return {'train': train, 'val': val, 'test': test}

    def compute_cluster_split(self, type, train_ratio, val_ratio, test_ratio, random_state):
        cluster_ids = [p['protein'][f'{type}_cluster_0.7'] for p in self.proteins]
        num_clusters = max(cluster_ids)
        test_size = int(num_clusters*test_ratio)
        val_size = int(num_clusters*val_ratio)
        unique, counts = np.unique(cluster_ids, return_counts=True)
        cluster_sizes = dict(zip(unique, counts))
        seq_threshold = int(np.median(list(cluster_sizes.values())))
        pool = [cluster for cluster, count in cluster_sizes.items() if count >= seq_threshold]
        np.random.seed(random_state)
        np.random.shuffle(pool)
        test_clusters, val_clusters = pool[:test_size], pool[test_size:test_size+val_size]
        inds = np.arange(self.size)
        def sample(c):
            return np.random.choice(inds[cluster_ids==c], size=seq_threshold, replace=False).tolist()
        test = list(itertools.chain.from_iterable(sample(c) for c in test_clusters))
        val = list(itertools.chain.from_iterable(sample(c) for c in val_clusters))
        train = [c for c in cluster_ids if c not in test and c not in val]
        return {'train': train, 'val': val, 'test': test}


    def compute_token_map(self):
        """ Computes and sets a dictionary that maps discrete labels to integers for classification tasks."""
        return None

    @property
    def task_type(self):
        """ Returns a string describing the type of task."""
        raise NotImplementedError

    @property
    def num_features(self):
        """ Number of input features to use for this task """
        raise NotImplementedError

    @property
    def num_classes(self):
        """ Size of the output dimension for this task """
        raise NotImplementedError

    @property
    def target(self, protein):
        """ Return the prediction target for one protein in the dataset.

        Parameters:
        ------------
        protein: dict
            proteinshake protein dictionary
        """
        raise NotImplementedError

    def evaluate(self, y_pred):
        """ Evaluates prediction quality.

        Parameters:
        -----------
        y_pred: list
            List of predicted outputs.

        Returns:
        --------
        dict
            Dictionary with evaluation results. Key-value pairs correspond to metric-score pairs. E.g. 'roc-auc': 0.7
        """
        raise NotImplementedError

    @property
    def train(self):
        return self.dataset[self.train_ind]

    @property
    def val(self):
        return self.dataset[self.val_ind]

    @property
    def test(self):
        return self.dataset[self.test_ind]

    def to_graph(self, *args, **kwargs):
        self.dataset = self.dataset.to_graph(*args, **kwargs)
        return self

    def to_point(self, *args, **kwargs):
        self.dataset = self.dataset.to_point(*args, **kwargs)
        return self

    def to_voxel(self, *args, **kwargs):
        self.dataset = self.dataset.to_voxel(*args, **kwargs)
        return self

    def pyg(self, *args, **kwargs):
        self.dataset = self.dataset.pyg(*args, **kwargs)
        return self

    def dgl(self, *args, **kwargs):
        self.dataset = self.dataset.dgl(*args, **kwargs)
        return self

    def nx(self, *args, **kwargs):
        self.dataset = self.dataset.nx(*args, **kwargs)
        return self

    def np(self, *args, **kwargs):
        self.dataset = self.dataset.np(*args, **kwargs)
        return self

    def tf(self, *args, **kwargs):
        self.dataset = self.dataset.tf(*args, **kwargs)
        return self

    def torch(self, *args, **kwargs):
        self.dataset = self.dataset.torch(*args, **kwargs)
        return self
