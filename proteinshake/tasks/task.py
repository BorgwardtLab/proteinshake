import os
import requests

import numpy as np
from sklearn.model_selection import train_test_split

from proteinshake.utils import save, load

class Task:
    """ Base class for a task.

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
    split: str, default='random'
        How to split the data. Can be 'random', 'sequence', 'structure' or 'custom'.
    split_similarity_threshold: float
        Maximum similarity between train and test set.
    """

    DatasetClass = None

    type = None
    input = None
    output = None
    default_metric = None

    def __init__(self,
                 root                       = 'data',
                 split                      = 'random',
                 split_similarity_threshold = 0.7,
                 use_precomputed_task       = True,
                 **kwargs
                ):
        self.root = root
        self.dataset = self.DatasetClass(root=root, **kwargs)
        self.split = split
        self.split_similarity_threshold = split_similarity_threshold

        if not self.files_exist:
            if use_precomputed_task and self.files_hosted: self.download_precomputed()
            else: self.compute()
        
        task_info = load(f'{self.root}/{self.__class__.__name__}.json.gz')
        split_info = task_info[f'{split}_split'][f'similarity {split_similarity_threshold}']

        self.token_map = task_info['token_map']
        self.targets = np.array(task_info['targets'])
        self.train_indices = np.array(split_info['train'])
        self.test_indices = np.array(split_info['test'])
        self.val_indices = np.array(split_info['val'])
        self.train_targets = np.array(self.targets[self.train_indices])
        self.test_targets = np.array(self.targets[self.test_indices])
        self.val_targets = np.array(self.targets[self.val_indices])

    def __getattr__(self, key):
        """ Captures method calls and forwards them to the dataset if they are a representation or framework conversion.
        """
        if key in ['to_graph','to_point','to_voxel','np','dgl','pyg','torch','nx','tf']:
            def proxy(*args, **kwargs):
                self.dataset = getattr(self.dataset, key)(*args, **kwargs)
                return self
            return proxy
        else:
            return getattr(self, key, lambda: None)
        
    @property
    def files_exist(self):
        return os.path.exists(f'{self.root}/{self.__class__.__name__}.json.gz')

    @property
    def files_hosted(self):
        return requests.head(f'{self.dataset.repository_url}/{self.__class__.__name__}.npz', timeout=5).status_code == 200
    
    def download_precomputed(self):
        raise NotImplementedError

    def compute(self):
        """ Computes the task file with splits, token maps, and targets.
        """
        self.token_map = self.compute_token_map()
        save({
            'random_split': self.compute_random_split(),
            'sequence_split': self.compute_sequence_split(),
            'structure_split': self.compute_structure_split(),
            'custom_split': self.compute_custom_split(),
            'token_map': self.token_map,
            'targets': self.compute_targets(),
        }, f'{self.root}/{self.__class__.__name__}.json.gz')

    def compute_random_split(self, seed=0):
        indices = np.arange(len(self.dataset.proteins()))
        train_indices, test_val_indices = train_test_split(indices, test_size=0.2, random_state=seed)
        test_indices, val_indices = train_test_split(test_val_indices, test_size=0.5, random_state=seed)
        return { 'similarity 0.7': {
            'train': train_indices.tolist(),
            'test': test_indices.tolist(),
            'val': val_indices.tolist(),
        }}

    def compute_sequence_split(self, seed=0):
        return None

    def compute_structure_split(self, seed=0):
        return None

    def compute_custom_split(self, seed=0):
        return None

    def compute_token_map(self):
        return None

    def compute_targets(self):
        return [
            self.target(protein_dict)
            for protein_dict in self.dataset.proteins()
        ]

    @property
    def train(self):
        return self.dataset[self.train_index]

    @property
    def val(self):
        return self.dataset[self.val_index]

    @property
    def test(self):
        return self.dataset[self.test_index]

    def target(self, protein_dict):
        """ Return the prediction target for one protein in the dataset.
        Can be a scalar or a fixed size array.

        Arguments
        ------------
        protein_dict: dict
            A protein dictionary.

        Returns
        -------
        target: float
            The prediction target.

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
    
    def dummy(self):
        raise NotImplementedError
    
            
