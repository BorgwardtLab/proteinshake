import os
import os.path as osp
import json
from pathlib import Path
from joblib import Parallel, delayed

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.cluster import AgglomerativeClustering

from proteinshake.utils import tmalign_wrapper, cdhit_wrapper

class ShakeTask:
    """ Base class for task-related utilities.
    This class wraps a proteinshake dataset and exposes split indices,
    integer-coded labels for classification tasks, and an evaluator function.
    Users should use this class to build their own dataloaders for training
    and evaluating models.

    Sample usage (assuming you have a model and dataloader in the namespace):

     .. code-block:: python

        >>> from proteinshake.tasks import EnzymeCommissionTask
        >>> task = EnzymeCommissionTask(root='foo')
        >>> task.train_ind
        ... [1, 3, 4, 5, 9, ...]
        >>> out = model(dataloader[task.test_ind])
        >>> targets = [task.target(i) for i in task.test_ind]
        >>> task.evaluate(out, targets)
        ... {'roc_auc_score': 0.7}


    Parameters
    ----------
    dataset: pytorch.datasets.Dataset
        Dataset to use for this task.
    random_state: int, default=42
        Random seed for reproducible splitting.
    train_ratio: float, default=.75
        Fraction of dataset to use for training.
    val_ratio: float, default=.15
        Fraction of dataset to use for validation.
    test_ratio: float, default=.10
        Fraction of dataset to use for testing.
    cace_dir: str, default='.proteinshake'
        Directory where we store the result of computing splits and tokenizing.
    """
    def __init__(self,
                 dataset,
                 random_state=42,
                 train_ratio=0.75,
                 val_ratio=.15,
                 test_ratio=.1,
                ):
        self.dataset = dataset
        self.train_ratio = train_ratio
        self.validation_ratio = val_ratio
        self.test_ratio = test_ratio
        self.random_state = random_state
        self.cache_dir = osp.join(self.dataset.root, "tasks")

        proteins, size = list(self.dataset.proteins())
        self.proteins = list(proteins)
        self.size = size

        self._process()

    def _process(self):
        """ Skeleton for processing a task. Tries to load results from previous output
        and then computes splits and label set.
        """
        self.create_cache()
        cache_dict = self.load_cache()
        if cache_dict is None:
            print(">>> computing task info.")
            self.compute_splits()
            try:
                self.compute_token_map()
            except NotImplementedError:
                print(">>> No tokenizer implemented. Make sure this is a regression task.")
                self.token_map = None
                pass
            self.process()
            self.cache()
        else:
            for k, v in cache_dict.items():
                setattr(self, k, v)

    def process(self):
        """ Override for extra processing"""
        pass

    def compute_splits(self):
        """ Compute train/val/test splits and sets respective attributes as lists of indices..
        """
        print(f">>> computing splits")
        inds = range(self.size)
        train, test = train_test_split(inds, test_size=1 - self.train_ratio)
        val, test = train_test_split(test, test_size=self.test_ratio/(self.test_ratio + self.validation_ratio))

        self.train_ind = train
        self.val_ind = val
        self.test_ind = test

    @property
    def train(self):
        return self.dataset[self.train_ind]

    @property
    def test(self):
        return self.dataset[self.test_ind]

    @property
    def val(self):
        return self.dataset[self.val_ind]

    def compute_token_map(self):
        """ Computes and sets a dictionary that maps discrete labels to integers for classification tasks."""
        raise NotImplementedError

    def label_process(self, label):
        """ Returns a cleaned prediction label. """
        return label

    @property
    def target(self, protein):
        """ Return the prediction target for one protein in the dataset.

        Parameters:
        ------------
        protein: dict
            proteinshake protein dictionary
        """
        raise NotImplementedError

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

    def evaluate(self, pred, true):
        """ Evaluates prediction quality.

        Parameters:
        -----------
        pred: list
            List of predicted outputs.
        true: list
            List of target values (output from the `ShakeTask.target()` method.

        Returns:
        --------
        dict
            Dictionary with evaluation results. Key-value pairs correspond to metric-score pairs. E.g. 'roc-auc': 0.7
        """
        raise NotImplementedError

    def create_cache(self):
        """ Creates the task info cache directory """
        if not os.path.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

    def load_cache(self):
        """ Tries to load the cache. Returns None if it does not exist

        Returns:
            dict:
                The dictionary storing task attributes. Returns None if no cache exists.
        """
        try:
            with open(osp.join(self.cache_dir, f"{self.dataset.name}.json"), "r") as t:
                task_dict = json.load(t)
                return task_dict
        except FileNotFoundError:
            return None

    def cache(self):
        """ Cache the stuff that needs to iterate over the dataset."""

        with open(osp.join(self.cache_dir, f"{self.dataset.name}.json"), "w") as t:
            json.dump({'train_ind': self.train_ind,
                       'test_ind': self.test_ind,
                       'val_ind': self.val_ind,
                       'token_map': self.token_map
                       }, t)
                       
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
