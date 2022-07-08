import os, shutil, argparse
from torch_pdb.datasets import PDBBindRefined, PDBBindPPI, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset

parser = argparse.ArgumentParser(description='Script to generate all datasets for release.')
parser.add_argument('--path', type=str, help='Path to store the final dataset objects.', default='.')
parser.add_argument('--scratch', type=str, help='Path to scratch (if on cluster).', default=os.path.expandvars('$SCRATCH'))
parser.add_argument('--njobs', type=int, help='Number of jobs.', default=20)
args = parser.parse_args()

PATH = args.path
SCRATCH = args.scratch if args.scratch != '' else args.path
n_jobs = args.njobs

os.makedirs(f'{PATH}/release', exist_ok=True)
for Dataset in [RCSBDataset, ECDataset, PDBBindRefined, PDBBindPPI,  TMScoreBenchmark, GODataset, PfamDataset]:
    name = Dataset.__name__
    print()
    print(name)
    ds = Dataset(root=f'{SCRATCH}/release/{name}', use_precomputed=False, n_jobs=n_jobs)
    del ds
    if SCRATCH != PATH and not os.path.exists(f'{PATH}/release/{name}.json.gz'):
        print('Copying...')
        shutil.copyfile(f'{SCRATCH}/release/{name}/{name}.json.gz', f'{PATH}/release/{name}.json.gz')

if SCRATCH != PATH and not os.path.exists(f'{PATH}/release/tmalign.json.gz'):
    print('Copying TM scores...')
    shutil.copyfile(f'{SCRATCH}/release/TMScoreBenchmark/tmalign.json.gz', f'{PATH}/release/tmalign.json.gz')


from torch_pdb.datasets.alphafold import AF_DATASET_NAMES
from torch_pdb.datasets import AlphaFoldDataset

for organism in AF_DATASET_NAMES.keys():
    print()
    print('AlphaFoldDataset', organism)
    ds = AlphaFoldDataset(root=f'{SCRATCH}/release/AlphaFoldDataset_{organism}', organism=organism, use_precomputed=False, n_jobs=n_jobs)
    del ds
    if SCRATCH != PATH and not os.path.exists(f'{PATH}/release/AlphaFoldDataset_{organism}.json.gz'):
        print('Copying...')
        shutil.copyfile(f'{SCRATCH}/release/AlphaFoldDataset_{organism}/AlphaFoldDataset.json.gz', f'{PATH}/release/AlphaFoldDataset_{organism}.json.gz')
