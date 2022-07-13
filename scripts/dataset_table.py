""" Print data statistics to stdout as markdown and latex tables
"""
import tempfile

import pandas as pd

from torch_pdb.datasets import PDBBindRefined, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset, AlphaFoldDataset
from torch_pdb.datasets.alphafold import AF_DATASET_NAMES

datasets = [
            RCSBDataset,
            PfamDataset,
            GODataset,
            ECDataset,
            PDBBindRefined,
            TMScoreBenchmark,
            ]

rows = []
for i, dataset in enumerate(datasets):
    with tempfile.TemporaryDirectory() as tmp:
        ds = dataset(root=tmp)
        desc = ds.describe()
        rows.append(desc)

# do alphafold seaprately
for org in AF_DATASET_NAMES.keys():
    print(org)
    with tempfile.TemporaryDirectory() as tmp:
        af_all = AlphaFoldDataset(root=tmp, organism=org)
        desc = af_all.describe()
        desc['name'] += f'_{org}'
        print(desc)
        rows.append(desc)
    df = pd.DataFrame(rows)
    md = df.to_markdown(index=False)
    tx = df.to_latex(index=False, na_rep='-')
    print(md)
    print()
    print(tx)


df = pd.DataFrame(rows)
md = df.to_markdown(index=False)
tx = df.to_latex(index=False)

print(md)
print()
print(tx)
