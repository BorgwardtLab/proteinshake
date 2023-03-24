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

    Sample usage (assuming you have a model in the namespace):

     .. code-block:: python

        >>> from proteinshake.tasks import EnzymeClassTask
        >>> task = EnzymeClassTask()
        >>> data = task.dataset.to_graph(eps=8).pyg()
        >>> y_pred = model(data[task.train])
        >>> task.evaluate(y_pred)
        ... {'roc_auc_score': 0.7}


    Arguments
    ----------
    dataset: pytorch.datasets.Dataset
        Dataset to use for this task.
    split: str, default='random'
        How to split the data. Can be 'random', 'sequence', or 'structure'.
    split_similarity_threshold: float
        Maximum similarity to allow between train and test samples.
    use_precomputed: bool, default=True
    """

    DatasetClass = None

    def __init__(self,
                 root                       = 'data',
                 split                      = 'random',
                 split_similarity_threshold = 0.7,
                 use_precomputed            = True,
                 **kwargs
                ):
        self.root = root
        self.dataset = self.DatasetClass(root=root, use_precomputed=use_precomputed)
        proteins = self.dataset.proteins()
        self.size = len(proteins)
        self.split_similarity_threshold = split_similarity_threshold
        self.split = split
    
        class Proteins(): # dummy class to implement __getitem__, could be implemented directly on the task
            def __init__(self, proteins):
                self.proteins = list(proteins)

            def __len__(self):
                return len(self.proteins)

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
        if not self.split == 'none':
            self.compute_index()
            self.compute_targets()

    def compute_index(self):
        split_name = f'{self.split}_split_{self.split_similarity_threshold}' if self.split in ['sequence','structure'] else f'{self.split}_split'
        if split_name in self.proteins[0]['protein']:
            self.train_index = np.array([i for i,p in enumerate(self.proteins) if p['protein'][split_name] == 'train'])
            self.val_index = np.array([i for i,p in enumerate(self.proteins) if p['protein'][split_name] == 'val'])
            self.test_index = np.array([i for i,p in enumerate(self.proteins) if p['protein'][split_name] == 'test'])
        else:
            self.train_index, self.val_index, self.test_index = self.compute_custom_split(self.split)

        self.update_index()

    def update_index(self):
        pass

    def compute_targets(self):
        # compute targets (e.g. for scaling)
        self.train_targets = np.array([self.target(self.proteins[i]) for i in self.train_index], dtype=object)
        self.val_targets = np.array([self.target(self.proteins[i]) for i in self.val_index], dtype=object)
        self.test_targets = np.array([self.target(self.proteins[i]) for i in self.test_index], dtype=object)
            
    def compute_custom_split(self, split):
        """ Implements custom splitting. Only necessary when not using the precomputed splits, e.g. when implementing a custom task.
        Note that the random, sequence and structure splits will be automatically computed for your custom task if it is merged into ProteinShake main.
        Compare also the proteinshake_release repository.

        Arguments
        ------------
        split: str
            Name of the custom split as passed to the task. ('random', 'sequence', 'structure', 'none')

        Returns:
        --------
        train_index
            Numpy array with the index of proteins in the train split.
        val_index
            Numpy array with the index of proteins in the validation split.
        test_index
            Numpy array with the index of proteins in the test split.
        """
        raise Exception('The requested split is not available. Implement <compute_custom_split> to provide your own splitting logic.')

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

        Arguments
        ------------
        protein: dict
            proteinshake protein dictionary


        .. code-block: python

            >>> from proteinshake.tasks import EnzymeCommissionTask
            >>> ta = EnzymeCommissionTask()
        """
        raise NotImplementedError

    def evaluate(self, y_true, y_pred):
        """ Evaluates prediction quality.

        Arguments
        -----------
        y_true: list
            List of ground truth outputs, (e.g. task.test_targets).
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
