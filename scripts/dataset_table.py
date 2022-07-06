""" Print data statistics to stdout as markdown and latex tables
"""
import tempfile

import pandas as pd

from torch_pdb.datasets import PDBBindRefined, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset


datasets = [
            RCSBDataset,
            PfamDataset,
            GODataset,
            ECDataset,
            PDBBindRefined,
            TMScoreBenchmark]

rows = []
for i, dataset in enumerate(datasets):
    ds = dataset(root=f"data_{i}", name='test')
    desc = ds.describe()
    print(desc)
    rows.append(desc)

df = pd.DataFrame(rows)
md = df.to_markdown(index=False)
tx = df.to_latex(index=False)

print(md)
print()
print(tx)
