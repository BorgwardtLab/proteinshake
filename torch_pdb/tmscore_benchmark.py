# -*- coding: utf-8 -*-
import glob, torch, requests, os, itertools, re, subprocess
from torch_geometric.data import extract_tar, download_url
from torch_pdb import TorchPDBDataset
from collections import defaultdict
from joblib import Parallel, delayed
from tqdm import tqdm

# short-term absolute path hack for TMalign
# we need to include this with the setuptools
import os
TMPATH = os.path.dirname(os.path.realpath(__file__))+'/../TMalign'

def tmalign_wrapper(pdb1, pdb2):
    try:
        out = subprocess.run([TMPATH,'-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return -1.
    return float(TM1), float(TM2), float(RMSD)


class TMScoreBenchmark(TorchPDBDataset):

    def __init__(self, use_precomputed=True, **kwargs):
        self.use_precomputed = use_precomputed
        super().__init__(**kwargs)
        self.tm_score, self.rmsd = self.compute_distances()

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')

    def get_id_from_filename(self, filename):
        return filename[:-4]

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
        if self.use_precomputed:
            download_url('https://github.com/BorgwardtLab/torch-pdb/releases/download/v1.0.0/tm-bench.tar.gz', f'{self.root}/raw')
            extract_tar(f'{self.root}/raw/tm-bench.tar.gz', f'{self.root}/raw')
        self.download_complete()

    def compute_distances(self):
        if os.path.exists(f'{self.root}/raw/tm-bench.pt'):
            return torch.load(f'{self.root}/raw/tm-bench.pt')

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

        torch.save((dist, rmsd), f'{self.root}/raw/tm-bench.pt')
        return dist, rmsd
