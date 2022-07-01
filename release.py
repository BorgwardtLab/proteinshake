'''
Script to generate all datasets for release.
'''

import os, shutil
from torch_pdb import PDBBindRefined, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset

PATH = os.path.dirname(os.path.realpath(__file__))
SCRATCH = os.path.expandvars('$SCRATCH') # this is for the cluster, set to PATH if run on local machine
n_jobs = 20

print('RCSBDataset')
ds = RCSBDataset(root=SCRATCH+'/release/rcsb', name='rcsb', n_jobs=n_jobs)
del ds

print('ECDataset')
ds = ECDataset(root=SCRATCH+'/release/ec', name='ec', n_jobs=n_jobs)
del ds

print('GODataset')
ds = GODataset(root=SCRATCH+'/release/go', name='go', n_jobs=n_jobs)
del ds

print('PfamDataset')
ds = PfamDataset(root=SCRATCH+'/release/pfam', name='pfam', n_jobs=n_jobs)
del ds

print('PDBBindRefined')
ds = PDBBindRefined(root=SCRATCH+'/release/pdbbind', name='pdbbind', n_jobs=n_jobs)
del ds

print('TMScoreBenchmark')
ds = TMScoreBenchmark(root=SCRATCH+'/release/tmscore', name='tmscore', use_precomputed=False, n_jobs=n_jobs)
del ds

if not PATH == SCRATCH:
    os.makedirs(f'{PATH}/release', exist_ok=True)
    for name in ['ec','go','pfam','rcsb','pdbbind','tmscore']:
        shutil.copyfile(f'{SCRATCH}/release/{name}/processed/{name}.pt', f'{PATH}/release/{name}.pt')
