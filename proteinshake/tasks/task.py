import os
import json
import itertools

import numpy as np
from sklearn.model_selection import train_test_split

from proteinshake.utils import download_url, save, load

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
                 root                       = 'data',
                 split                      = 'random',
                 split_similarity_threshold = 0.7,
                 random_state               = 42,
                 train_ratio                = 0.80,
                 val_ratio                  = 0.10,
                 test_ratio                 = 0.10,
                 use_precomputed            = True,
                 **kwargs
                ):
        self.root = root
        self.dataset = self.DatasetClass(root=root)
        proteins = self.dataset.proteins()
        self.size = len(proteins)
        class Proteins(): # dummy class to implement __getitem__, could be implemented directly on the task
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
        self.proteins = Proteins(proteins)
        self.name = self.__class__.__name__

        # load split indices
        split_name = f'{split}_split_{split_similarity_threshold}' if split in ['sequence','structure'] else f'{split}_split'
        if split_name in self.proteins[0]['protein']:
            self.train_index = np.array([i for i,p in enumerate(self.proteins) if p['protein'][split_name] == 'train'])
            self.val_index = np.array([i for i,p in enumerate(self.proteins) if p['protein'][split_name] == 'test'])
            self.test_index = np.array([i for i,p in enumerate(self.proteins) if p['protein'][split_name] == 'val'])
        else:
            self.train_index, self.val_index, self.test_index = self.compute_custom_split(split)

        # compute targets (e.g. for scaling)
        self.train_targets = np.array([self.target(self.proteins[i]) for i in self.train_index])
        self.val_targets = np.array([self.target(self.proteins[i]) for i in self.val_index])
        self.test_targets = np.array([self.target(self.proteins[i]) for i in self.test_index])
            
    def compute_custom_split(self, split):
        """ Implements custom splitting. Only necessary when not using the precomputed splits, e.g. when implementing a custom task.
        Note that the random, sequence and structure splits will be automatically computed for your custom task if it is merged into ProteinShake main.
        Compare also the proteinshake_release repository.

        Parameters:
        ------------
        split: str
            Name of the custom split as passed to the task.

        Returns:
        --------
        train_index
            Numpy array with the index of proteins in the train split.
        val_index
            Numpy array with the index of proteins in the validation split.
        test_index
            Numpy array with the index of proteins in the test split.
        """

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
        return self.dataset[self.train_index]

    @property
    def val(self):
        return self.dataset[self.val_index]

    @property
    def test(self):
        return self.dataset[self.test_index]

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
