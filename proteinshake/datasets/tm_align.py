# -*- coding: utf-8 -*-
'''
TMalign needs to be in your $PATH. Follow the instructions at https://zhanggroup.org/TM-align/readme.c++.txt
'''
import glob
import requests
import os
import itertools
import re
import subprocess
from collections import defaultdict

from joblib import Parallel, delayed
from tqdm import tqdm

from proteinshake.datasets import Dataset
from proteinshake.utils import (extract_tar,
                                download_url,
                                save,
                                load,
                                unzip_file,
                                tmalign_wrapper
                                )


class TMAlignDataset(Dataset):
    """Proteins with TM scores between all pairs.
    This dataset contains 200 proteins from the TMalign benchmark dataset.
    The dataset has a global attribute `tm_score` which is a dictionary
    containing the TMscore between all pairs of proteins in the dataset.

    .. code-block:: python

        from proteinshake.datasets import TMScoreBenchmark

        dataset = TMScoreBenchmark()
        protein_1, protein_2 = dataset[0].name, dataset[2].name

        dataset[protein_1][protein_2]
        >>> 0.32

    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.
    use_precomputed: bool
        If `True` uses TM scores from saved TMalign output. Otherwise, recomputes.
    """

    def get_raw_files(self, compressed=False):
        return glob.glob(f'{self.root}/raw/files/*.pdb')

    def get_id_from_filename(self, filename):
        return filename[:-4]

    def download(self):
        lines = requests.get("https://zhanggroup.org/TM-align/benchmark/").text.split("\n")
        links = []
        print('Downloading TMScore Benchmark PDBs...')
        for l in lines:
            m = re.search(".pdb", l)
            if m:
                start, end = m.span()
                pdbid = l[start-5:end]
                links.append(f"https://zhanggroup.org/TM-align/benchmark/{pdbid}")
        links = links[:self.limit] # for testing
        for link in tqdm(links):
            download_url(link, f'{self.root}/raw/files', log=False)

    def describe(self):
        desc = super().describe()
        desc['property'] = 'TM Score'
        desc['values'] = "[0-1]"
        desc['type'] = "Real-valued, Pairwise"
        return desc
