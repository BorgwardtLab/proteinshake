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

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.tm_score, self.rmsd = self.compute_distances()

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
        links = links[:self.download_limit()] # for testing
        for link in tqdm(links):
            download_url(link, f'{self.root}/raw/files', log=False)

    def compute_distances(self, n_jobs=1):
        """ Launch TMalign on all pairs of proteins in dataset.
        Saves TMalign output to `self.raw_dir/tmalign.json.gz`
        Returns
        -------
        dict
            TM score between all pairs of proteins as a dictionary.
        dict
            RMSD between all pairs of proteins as a dictionary.
        """
        if os.path.exists(f'{self.root}/tmalign.json'):
            return load(f'{self.root}/tmalign.json')
        elif self.use_precomputed:
            download_url(f'{self.repository_url}/tmalign.json.gz', f'{self.root}')
            print('Unzipping...')
            unzip_file(f'{self.root}/tmalign.json.gz')
            return load(f'{self.root}/tmalign.json')
        if self.n_jobs == 1:
            print('Computing the TM scores with use_precompute = False is very slow. Consider increasing n_jobs.')

        pdbs = self.get_raw_files()
        pairs = list(itertools.combinations(range(len(pdbs)), 2))
        todo = [(pdbs[p1], pdbs[p2]) for p1, p2 in pairs]

        output = Parallel(n_jobs=self.n_jobs)(
            delayed(tmalign_wrapper)(*pair) for pair in tqdm(todo, desc='Computing TM Scores')
        )

        dist = defaultdict(lambda: {})
        rmsd = defaultdict(lambda: {})
        for (pdb1, pdb2), d in zip(todo, output):
            name1 = os.path.basename(pdb1).split('.')[0]
            name2 = os.path.basename(pdb2).split('.')[0]
            dist[name1] = {**dist[name1], name2: d[0]}
            dist[name2] = {**dist[name2], name1: d[1]}
            rmsd[name1] = {**rmsd[name1], name2: d[2]}
        dist = dict(dist)
        rmsd = dict(rmsd)

        save((dist, rmsd), f'{self.root}/tmalign.json')
        return dist, rmsd

    def describe(self):
        desc = super().describe()
        desc['property'] = 'TM Score'
        desc['values'] = "[0-1]"
        desc['type'] = "Real-valued, Pairwise"
        return desc
