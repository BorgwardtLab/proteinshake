import os, shutil, argparse
from datetime import datetime
from proteinshake.utils import zip_file
import importlib
from tqdm import tqdm
from functools import partialmethod

# disable tqdm
tqdm.__init__ = partialmethod(tqdm.__init__, disable=True)

datasets = importlib.import_module('proteinshake.datasets')

parser = argparse.ArgumentParser(description='Script to generate all datasets for release.')
parser.add_argument('--path', type=str, help='Path to store the final dataset objects.', default='./release')
parser.add_argument('--scratch', type=str, help='Path to scratch (if on cluster).', default=os.path.expandvars(f'$TMPDIR/proteinshake'))
parser.add_argument('--tag', type=str, help='Release tag', default=datetime.now().strftime('%d%b%Y').upper())
parser.add_argument('--njobs', type=int, help='Number of jobs.', default=10)
parser.add_argument('--dataset', type=str, help='Name of the dataset class (case sensitive)', default='RCSBDataset')
parser.add_argument('--organism', type=str, help='Organism (for AlphaFold datasets)', default='swissprot')
args = parser.parse_args()

RELEASE = args.tag
PATH = args.path+'/'+RELEASE
SCRATCH = args.scratch+'/'+RELEASE if args.scratch != '' else args.path
DATASET = args.dataset
ORGANISM = args.organism
n_jobs = args.njobs
clustering = not DATASET in ['RCSBDataset','AlphaFoldDataset']

if os.path.exists(f'{PATH}/{DATASET}.atom.avro.gz'):
    print(f'Skipping {DATASET}')
    exit()

print()
print(DATASET)

os.makedirs(f'{SCRATCH}', exist_ok=True)
os.makedirs(f'{PATH}', exist_ok=True)

# create dataset
Dataset = getattr(datasets, DATASET)
kwargs = {}
if DATASET == 'AlphaFoldDataset':
    kwargs['organism'] = ORGANISM
    DATASET = f'{DATASET}_{ORGANISM}'
kwargs = {
    'root'                            : f'{SCRATCH}/{DATASET}',
    'use_precomputed'                 : False,
    'n_jobs'                          : n_jobs,
    'cluster_structure'               : clustering,
    'cluster_sequence'                : clustering,
    'similarity_threshold_structure'  : [0.5, 0.6, 0.7, 0.8, 0.9],
    'similarity_threshold_sequence'   : [0.5, 0.6, 0.7, 0.8, 0.9],
    **kwargs
}
ds = Dataset(**kwargs)

# zipping
print('Compressing...')
zip_file(f'{SCRATCH}/{DATASET}/{DATASET}.residue.avro')
zip_file(f'{SCRATCH}/{DATASET}/{DATASET}.atom.avro')
if clustering:
    zip_file(f'{SCRATCH}/{DATASET}/{DATASET}.cdhit.json')
    zip_file(f'{SCRATCH}/{DATASET}/{DATASET}.tmalign.json')

# copy from scratch to target
print('Copying...')
shutil.copyfile(f'{SCRATCH}/{DATASET}/{DATASET}.residue.avro.gz', f'{PATH}/{DATASET}.residue.avro.gz')
shutil.copyfile(f'{SCRATCH}/{DATASET}/{DATASET}.atom.avro.gz', f'{PATH}/{DATASET}.atom.avro.gz')
if clustering:
    shutil.copyfile(f'{SCRATCH}/{DATASET}/{DATASET}.cdhit.json.gz', f'{PATH}/{DATASET}.cdhit.json.gz')
    shutil.copyfile(f'{SCRATCH}/{DATASET}/{DATASET}.tmalign.json.gz', f'{PATH}/{DATASET}.tmalign.json.gz')

# cleaning
print('Cleaning up...')
shutil.rmtree(f'{SCRATCH}/{DATASET}')
