# -*- coding: utf-8 -*-
import glob
import torch
import requests
import os
import itertools
import re
import subprocess
from collections import defaultdict
from joblib import Parallel, delayed

from tqdm import tqdm
from torch_geometric.data import extract_tar, download_url

from torch_pdb.datasets import TorchPDBDataset

# short-term absolute path hack for TMalign
# we need to include this with the setuptools
TMPATH = os.path.dirname(os.path.realpath(__file__))+'/../TMalign'

def tmalign_wrapper(pdb1, pdb2):
    """Compute TM score with TMalign between two PDB structures.

    Parameters
    ----------
    pdb1: str
        Path to PDB.
    arg2 : str
        Path to PDB.

    Returns
    -------
    float
        TM score from `pdb1` to `pdb2`
    float
        TM score from `pdb2` to `pdb1`
    float
        RMSD between structures
    """
    try:
        out = subprocess.run([TMPATH,'-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return -1.
    return float(TM1), float(TM2), float(RMSD)


class TMScoreBenchmark(TorchPDBDataset):
    """"Dataset conatining TM scores between pairs of proteins.

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

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')

    def get_id_from_filename(self, filename):
        return filename[:-4]

    def download_precomputed(self):
        super().download_precomputed()
        if self.use_precomputed:
            download_url(f'https://github.com/BorgwardtLab/torch-pdb/releases/download/{self.release}/tm-bench.pt', f'{self.root}')

    def download(self):
        lines = requests.get("https://zhanggroup.org/TM-align/benchmark/").text.split("\n")
        pdblist = []
        print('Downloading TMScore Benchmark PDBs...')
        for l in lines:
            m = re.search(".pdb", l)
            if m:
                start, end = m.span()
                pdbid = l[start-5:end]
                download_url(f"https://zhanggroup.org/TM-align/benchmark/{pdbid}", f'{self.root}/raw/files', log=False)

    def compute_distances(self, n_jobs=1):
        """ Launch TMalign on all pairs of proteins in dataset.
        Saves TMalign output to `self.raw_dir/tm-bench.pt`

        Returns
        -------
        dict
            TM score between all pairs of proteins as a dictionary.
        dict
            RMSD between all pairs of proteins as a dictionary.
        """
        if os.path.exists(f'{self.root}/tm-bench.pt'):
            return torch.load(f'{self.root}/tm-bench.pt')
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

        torch.save((dist, rmsd), f'{self.root}/tm-bench.pt')
        return dist, rmsd