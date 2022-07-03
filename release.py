import os, shutil, argparse
from torch_pdb import PDBBindRefined, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset

parser = argparse.ArgumentParser(description='Script to generate all datasets for release.')
parser.add_argument('--path', type=str, help='Path to store the final dataset objects.', default='.')
parser.add_argument('--scratch', type=str, help='Path to scratch (if on cluster).', default=os.path.expandvars('$SCRATCH'))
parser.add_argument('--njobs', type=int, help='Number of jobs.', default=20)
args = parser.parse_args()

PATH = args.path
SCRATCH = args.scratch if args.scratch != '' else args.path
n_jobs = args.njobs

os.makedirs(f'{PATH}/release', exist_ok=True)
for Dataset in [RCSBDataset, ECDataset, PDBBindRefined, TMScoreBenchmark, GODataset, PfamDataset]:
    name = Dataset.__name__
    print()
    print(name)
    ds = Dataset(root=f'{SCRATCH}/release/{name}', name=name, use_precomputed=False, n_jobs=n_jobs)
    del ds
    if SCRATCH != PATH and not os.path.exists(f'{PATH}/release/{name}.pt'):
        print('Copying...')
        shutil.copyfile(f'{SCRATCH}/release/{name}/processed/{name}.pt', f'{PATH}/release/{name}.pt')
