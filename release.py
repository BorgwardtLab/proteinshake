import os
from torch_pdb import PDBBindRefined, TMScoreBenchmark, GODataset, ECDataset, PfamDataset, RCSBDataset

PATH = os.path.dirname(os.path.realpath(__file__))

print('RCSBDataset')
ds = RCSBDataset(root=PATH+'/data/rcsb', name='test')
print('test:', ds[0])
print(len(ds))
print()

print('PfamDataset')
ds = PfamDataset(root=PATH+'/data/pfam', name='test')
print('test:', ds[0]['Pfam'])
print(len(ds))
print()

print('ECDataset')
ds = ECDataset(root=PATH+'/data/ec', name='test')
print('test:', ds[0]['EC'])
print(len(ds))
print()

print('GODataset')
ds = GODataset(root=PATH+'/data/go', name='test')
print('test:', ds[0]['GO'])
print(len(ds))
print()

print('PDBBindRefined')
ds = PDBBindRefined(root=PATH+'/data/pdbbind', name='test')
print('test:', ds[0])
print(len(ds))
print()

print('TMScoreBenchmark')
ds = TMScoreBenchmark(root=PATH+'/data/tmscore', name='test', use_precomputed=False)
print('test:', ds.tm_score['1adr_']['1aho_'])
print(len(ds))
print()
