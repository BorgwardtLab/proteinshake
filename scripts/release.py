'''
Script to generate a dataset release for hosting. It runs all available datasets with 'use_precomputed=False'. If a scratch location (--scratch) is passed (default), the raw dataset will be downloaded to the scratch and afterwards copied to the target folder (--path).

How to generate a new release:
1. run release.py
2. make a new release on GitHub, attach the .json.gz file of each dataset
3. tag it after date of download, like '12JUL2022'
4. change the default argument 'release' in proteinshake.datasets.dataset.Dataset to the new release tag

How to add a new dataset to the release pipeline:
1. import the dataset class
2. add it to the first loop
3. run release

'''

import os, shutil, argparse
from datetime import datetime
from proteinshake.datasets import RCSBDataset, GeneOntologyDataset, EnzymeCommissionDataset, PfamDataset, ProteinProteinInterfaceDataset, ProteinLigandInterfaceDataset, TMAlignDataset, SCOPDataset
from proteinshake.utils import zip_file

parser = argparse.ArgumentParser(description='Script to generate all datasets for release.')
parser.add_argument('--path', type=str, help='Path to store the final dataset objects.', default='./release')
parser.add_argument('--scratch', type=str, help='Path to scratch (if on cluster).', default=os.path.expandvars(f'$SCRATCH/proteinshake'))
parser.add_argument('--tag', type=str, help='Release tag', default=datetime.now().strftime('%d%b%Y').upper())
parser.add_argument('--njobs', type=int, help='Number of jobs.', default=20)
args = parser.parse_args()

RELEASE = args.tag
PATH = args.path+'/'+RELEASE
SCRATCH = args.scratch+'/'+RELEASE if args.scratch != '' else args.path
n_jobs = args.njobs

os.makedirs(f'{SCRATCH}', exist_ok=True)
os.makedirs(f'{PATH}', exist_ok=True)

###################
# PDB Datasets
###################
for Dataset in [RCSBDataset, GeneOntologyDataset, EnzymeCommissionDataset, PfamDataset, ProteinProteinInterfaceDataset, ProteinLigandInterfaceDataset, TMAlignDataset, SCOPDataset]:
    # name it after class name
    name = Dataset.__name__
    if os.path.exists(f'{SCRATCH}/{name}.atom.avro.gz'):
        print(f'Skipping {name}')
        continue
    print()
    print(name)
    # create dataset
    ds = Dataset(root=f'{SCRATCH}/{name}', use_precomputed=False, n_jobs=n_jobs)
    print('Compressing...')
    zip_file(f'{SCRATCH}/{name}/{name}.residue.avro')
    zip_file(f'{SCRATCH}/{name}/{name}.atom.avro')
    # delete to free memory
    del ds
    # copy from scratch to target
    print('Copying...')
    shutil.copyfile(f'{SCRATCH}/{name}/{name}.residue.avro.gz', f'{PATH}/{name}.residue.avro.gz')
    shutil.copyfile(f'{SCRATCH}/{name}/{name}.atom.avro.gz', f'{PATH}/{name}.atom.avro.gz')
    #if SCRATCH != PATH and not os.path.exists(f'{PATH}/{name}.json.gz'):
    #    print('Copying...')
    #    shutil.copyfile(f'{SCRATCH}/{name}/{name}.json.gz', f'{PATH}/{name}.json.gz')
    if name == 'TMAlignDataset':
        # copy the extra file (pairwise distances) from the TM dataset
        print('Copying TM scores...')
        zip_file(f'{SCRATCH}/TMAlignDataset/tmalign.json')
        shutil.copyfile(f'{SCRATCH}/TMAlignDataset/tmalign.json.gz', f'{PATH}/tmalign.json.gz')
    # cleanup
    print('Cleaning up...')
    shutil.rmtree(f'{SCRATCH}/{name}')

###################
# AlphaFold Datasets
###################
from proteinshake.datasets.alphafold import AF_DATASET_NAMES
from proteinshake.datasets import AlphaFoldDataset

# create one dataset for each organism
for organism in AF_DATASET_NAMES.keys():
    if os.path.exists(f'{SCRATCH}/AlphaFoldDataset_{organism}.atom.avro.gz'):
        print(f'Skipping AlphaFold {organism}')
        continue
    print()
    print('AlphaFoldDataset', organism)
    ds = AlphaFoldDataset(root=f'{SCRATCH}/AlphaFoldDataset_{organism}', organism=organism, use_precomputed=False, n_jobs=n_jobs)
    print('Compressing...')
    zip_file(f'{SCRATCH}/AlphaFoldDataset_{organism}/AlphaFoldDataset.residue.avro')
    zip_file(f'{SCRATCH}/AlphaFoldDataset_{organism}/AlphaFoldDataset.atom.avro')
    del ds
    print('Copying...')
    shutil.copyfile(f'{SCRATCH}/AlphaFoldDataset_{organism}/AlphaFoldDataset.residue.avro.gz', f'{PATH}/AlphaFoldDataset_{organism}.residue.avro.gz')
    shutil.copyfile(f'{SCRATCH}/AlphaFoldDataset_{organism}/AlphaFoldDataset.atom.avro.gz', f'{PATH}/AlphaFoldDataset_{organism}.atom.avro.gz')
    # cleanup
    print('Cleaning up...')
    shutil.rmtree(f'{SCRATCH}/AlphaFoldDataset_{organism}')
    #if SCRATCH != PATH and not os.path.exists(f'{PATH}/AlphaFoldDataset_{organism}.json.gz'):
    #    print('Copying...')
    #    shutil.copyfile(f'{SCRATCH}/AlphaFoldDataset_{organism}/AlphaFoldDataset.json.gz', f'{PATH}/AlphaFoldDataset_{organism}.json.gz')
