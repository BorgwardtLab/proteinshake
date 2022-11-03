# -*- coding: utf-8 -*-

import os

import numpy as np
from tqdm import tqdm
from scipy import sparse

from proteinshake.utils import load, save
from proteinshake.frameworks.dataset import FrameworkDataset

class NumpyVoxelDataset(FrameworkDataset):
    def convert_to_framework(self, data_item):
        return sparse.csr_array(data_item.data)

    def load_transform(self, data, protein_dict):
        return data.todense(), protein_dict


class NumpyPointDataset(FrameworkDataset):
    def convert_to_framework(self, data_item):
        return data_item.data

    def load_transform(self, data, protein_dict):
        return data, protein_dict
