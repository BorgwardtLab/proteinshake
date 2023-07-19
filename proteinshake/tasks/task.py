import os
import requests
import itertools

import numpy as np
from sklearn.model_selection import train_test_split

from proteinshake.utils import save, load, progressbar

class Task:
    """ Base class for a task.

    Arguments
    ----------
    root: str, default 'data'
        The root directory to put the dataset and task data.
    split: str, default 'random'
        How to split the data. Can be 'random', 'sequence', 'structure' or 'custom'.
    split_similarity_threshold: float
        Maximum similarity between train and test set.
    """

    DatasetClass = None

    type = None
    input = None
    output = None
    default_metric = None
    pairwise = False
    level = 'Protein'

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
        
        task_data = load(f'{self.root}/{self.__class__.__name__}.json.gz')
        self.load_splits(task_data)
        self.load_targets(task_data)

    def __getattr__(self, key):
        """ Captures method calls and forwards them to the dataset if they are a representation or framework conversion.
        """
        if key in ['to_graph','to_point','to_voxel','np','dgl','pyg','torch','nx','tf']:
            def proxy(*args, **kwargs):
                self.dataset = getattr(self.dataset, key)(*args, **kwargs)
                return self
            return proxy
        else: return object.__getattribute__(self, key)
        
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
        self.targets = self.compute_paired_targets() if self.pairwise else self.compute_targets()
        save({
            'random_split': self.compute_random_split(),
            'sequence_split': self.compute_sequence_split(),
            'structure_split': self.compute_structure_split(),
            'custom_split': self.compute_custom_split(),
            'token_map': self.token_map,
            'targets': self.targets,
        }, f'{self.root}/{self.__class__.__name__}.json.gz')

    def compute_random_split(self):
        indices = np.arange(len(self.dataset.proteins()))
        train_indices, test_val_indices = train_test_split(indices, test_size=0.2, random_state=0)
        test_indices, val_indices = train_test_split(test_val_indices, test_size=0.5, random_state=0)
        return {
            'train': train_indices.tolist(),
            'test': test_indices.tolist(),
            'val': val_indices.tolist(),
        }

    def compute_sequence_split(self):
        return None

    def compute_structure_split(self):
        return None

    def compute_custom_split(self):
        return None

    def compute_token_map(self):
        return None

    def compute_targets(self):
        return [
            self.target(protein_dict)
            for protein_dict in progressbar(self.dataset.proteins(), verbosity=self.dataset.verbosity, desc='Computing targets')
        ]
        
    def compute_paired_targets(self):
        print(len(self.dataset.proteins()))
        index = self.compute_paired_index(np.arange(len(self.dataset.proteins())))
        print(len(index))
        exit()
        return [
            self.target(A,B)
            for A in progressbar(self.dataset.proteins(), verbosity=self.dataset.verbosity, desc='Computing targets')
            for B in self.dataset.proteins()
        ]
        
    def load_splits(self, task_data):
        split_data = task_data[f'{self.split}_split']
        if self.split != 'random': split_data = split_data[f'similarity_{self.split_similarity_threshold}']
        self.train_index, self.test_index, self.val_index = split_data['train'], split_data['test'], split_data['val']
        if self.pairwise:
            self.train_index = self.compute_paired_index(self.train_index)
            self.val_index = self.compute_paired_index(self.val_index)
            self.test_index = self.compute_paired_index(self.test_index)
    
    def compute_paired_index(self, index):
        """ Computes all pairs between each element in the index.
        """
        combinations = np.array(list(itertools.combinations(range(len(index)), 2)), dtype=int)
        return np.array(index)[combinations]
        
    def load_targets(self, task_data):
        self.token_map = task_data['token_map']
        self.targets = np.array([self.target_transform(target) for target in task_data['targets']], dtype=object)
        if self.pairwise: self.targets = self.targets.reshape(int(np.square(len(self.targets))), int(np.square(len(self.targets))))
        
    def target_transform(self, target):
        return target

    @property
    def X_train(self):
        return self.dataset[self.train_index]

    @property
    def X_test(self):
        return self.dataset[self.test_index]
    
    @property
    def X_val(self):
        return self.dataset[self.val_index]
    
    @property
    def y_train(self):
        return np.take(self.targets, self.train_index)

    @property
    def y_test(self):
        return np.take(self.targets, self.test_index)
    
    @property
    def y_val(self):
        return np.take(self.targets, self.val_index)
    
    @property
    def y_dummy(self):
        """ Generates a random data object of the right type and shape to compare against the test set.
        """
        raise NotImplementedError

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
    
    