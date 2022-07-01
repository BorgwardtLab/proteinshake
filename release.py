'''
Script to generate all datasets for release.
'''

import os, shutil
from torch_pdb import PDBBindRefined, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset

PATH = os.path.dirname(os.path.realpath(__file__))
n_jobs = 10

print('ECDataset')
ds = ECDataset(root=PATH+'/release/ec', name='ec', n_jobs=n_jobs)
del ds

print('GODataset')
ds = GODataset(root=PATH+'/release/go', name='go', n_jobs=n_jobs)
del ds

print('PfamDataset')
ds = PfamDataset(root=PATH+'/release/pfam', name='pfam', n_jobs=n_jobs)
del ds

print('RCSBDataset')
ds = RCSBDataset(root=PATH+'/release/rcsb', name='rcsb', n_jobs=n_jobs)
del ds

print('PDBBindRefined')
ds = PDBBindRefined(root=PATH+'/release/pdbbind', name='pdbbind', n_jobs=n_jobs)
del ds

print('TMScoreBenchmark')
ds = TMScoreBenchmark(root=PATH+'/release/tmscore', name='tmscore', use_precomputed=False, n_jobs=n_jobs)
del ds
