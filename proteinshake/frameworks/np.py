# -*- coding: utf-8 -*-

import os

import numpy as np
from tqdm import tqdm
from scipy import sparse

from proteinshake.utils import load, save
from proteinshake.frameworks.dataset import FrameworkDataset

class NumpyVoxelDataset(FrameworkDataset):
    """ Voxel dataset for NumPy.
    """

    def convert_to_framework(self, data_item):
        return data_item.data


class NumpyPointDataset(FrameworkDataset):
    """ Point dataset for NumPy.
    """

    def convert_to_framework(self, data_item):
        return data_item.data
