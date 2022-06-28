"""
Precompute TM scores on random pairs of graphs and dump as pickle.
Adds an attribute to the dataset object which is a dictionary containing
TM score for all pairs of data entries.
"""


import os
import tempfile
import gzip
import uuid
import shutil
import pickle
import random
import itertools
import multiprocessing as mp
import subprocess
import argparse
import requests
from collections import defaultdict

import wget
from tqdm import tqdm
from biopandas.pdb import PandasPdb
import numpy as np
# import seaborn as sns
import matplotlib.pyplot as plt
from joblib import Parallel, delayed

PATH = os.path.realpath(os.path.dirname(__file__))

def tmscore(pdb1, pdb2, n_jobs=1):
    """ Compute TM score for a pair of PDBs.
    Make sure `TMalign` is in `PATH`
    """
    try:
        out = subprocess.run(['TMalign','-outfmt','2', pdb1, pdb2], stdout=subprocess.PIPE).stdout.decode()
        path1, path2, TM1, TM2, RMSD, ID1, ID2, IDali, L1, L2, Lali = out.split('\n')[1].split('\t')
    except Exception as e:
        print(e)
        return -1.
    return float(TM1), float(TM2), float(RMSD)


def cline():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", type=str, help="Path to output folder, will be used for torch-pdb dataset construction.")
    parser.add_argument("--pdblist", type=str, help="Path to list of PDBs to use.")
    parser.add_argument("--custom-urls", type=str, default=None, help="Path to list of URLs to download. One line for each PDB..")
    parser.add_argument("--seed", type=int, default=42, help="random seed")
    parser.add_argument("--n-jobs", type=int, default=4, help="number of workers")

    args = parser.parse_known_args()[0]
    random.seed(args.seed)
    return args

def launch(args):
    """ Compute TM score for all pairs of proteins in the chosen dataset.
    """
    # fetch the PDBs
    pdb_dump = os.path.join(args.dest, 'pdbs')
    try:
        os.mkdir(args.dest)
        os.mkdir(pdb_dump)
    except FileExistsError:
        pass

    if not args.custom_urls is None:
        with open(args.custom_urls, 'r') as urls:
            for i, url in enumerate(urls):
                if i > 10:
                    break
                r = wget.download(url.strip(), out=pdb_dump)
    else:
        for pdbid in args.pdblist:
            ppdb = PandasPdb().fetch_pdb(pdbid)
            ppdb.to_pdb(os.path.join(pdb_dump, pdbid))

    pdbs = [os.path.join(pdb_dump, p) for p in os.listdir(pdb_dump)]
    pairs = list(itertools.combinations(range(len(pdbs)), 2))

    todo = [(pdbs[p1], pdbs[p2]) for p1, p2 in pairs]

    # compute distances
    distances = Parallel(n_jobs=args.n_jobs)(
        delayed(tmscore)(*pair) for pair in tqdm(todo, desc='TM Score')
    )
    dist_dict = defaultdict(lambda: {})
    rmsd = defaultdict(lambda: {})
    for (pdb1, pdb2), d in zip(todo, distances):
        name1 = os.path.basename(pdb1)
        name2 = os.path.basename(pdb2)
        dist_dict[name1] = {**dist_dict[name1], name2: d[0]}
        dist_dict[name2] = {**dist_dict[name2], name1: d[1]}
        rmsd[name1] = {**rmsd[name1], name2: d[2]}
    dist_dict = dict(dist_dict)
    rmsd = dict(rmsd)


    pickle.dump(dist_dict, open(os.path.join(args.dest, 'tm_scores.pkl'), "wb"))
    pickle.dump(rmsd, open(os.path.join(args.dest, "rmsd.pkl"), "wb"))

if __name__ == "__main__":
    args = cline()
    launch(args)
