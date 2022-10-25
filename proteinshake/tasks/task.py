import os
import os.path as osp
import json

from sklearn import metrics
from sklearn.model_selection import train_test_split


class ShakeTask:
    def __init__(self,
                 dataset,
                 random_state=42,
                 train_ratio=0.75,
                 val_ratio=.15,
                 test_ratio=.1,
                 cache_dir='.proteinshake'):
        self.dataset = dataset
        self.data_name = dataset.__class__.__name__
        self.train_ratio = train_ratio
        self.validation_ratio = val_ratio
        self.test_ratio = test_ratio
        self.random_state = random_state
        self.cache_dir = cache_dir

        proteins, size = list(self.dataset.proteins())
        self.proteins = list(proteins)
        self.size = size

        self._process()

    def _process(self):
        self.create_cache()
        cache_dict = self.load_cache()
        if cache_dict is None:
            print(">>> computing task info.")
            self.compute_splits()
            self.compute_token_map()
            self.process()
            self.cache()
        else:
            for k, v in cache_dict.items():
                setattr(self, k, v)

    def process(self):
        """ Override for extra processing"""
        pass

    def compute_splits(self):
        print(f">>> computing splits")
        inds = range(self.size)
        train, test = train_test_split(inds, test_size=1 - self.train_ratio)
        val, test = train_test_split(test, test_size=self.test_ratio/(self.test_ratio + self.validation_ratio))

        self.train_ind = train
        self.val_ind = val
        self.test_ind = test

    def compute_token_map(self):
        self.token_map = {}

    def label_process(self, label):
        """ Returns a cleaned prediction label. """
        return label

    @property
    def level(self):
        raise NotImplementedError

    @property
    def target(self, idx):
        raise NotImplementedError

    @property
    def level(self):
        raise NotImplementedError

    @property
    def num_features(self):
        raise NotImplementedError

    @property
    def num_labels(self):
        raise NotImplementedError

    def evaluator(self):
        raise NotImplementedError

    def create_cache(self):
        if not osp.exists(self.cache_dir):
            os.makedirs(self.cache_dir)

    def load_cache(self):
        try:
            with open(osp.join(self.cache_dir, f"{self.data_name}.json"), "r") as t:
                task_dict = json.load(t)
                return task_dict
        except FileNotFoundError:
            return None

    def cache(self):
        """ Cache the stuff that needs to iterate over the dataset."""

        with open(osp.join(self.cache_dir, f"{self.data_name}.json"), "w") as t:
            json.dump({'train_ind': self.train_ind,
                       'test_ind': self.test_ind,
                       'val_ind': self.val_ind,
                       'token_map': self.token_map
                       }, t)
