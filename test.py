import os
from torch_pdb import PDBBindRefined, TMScoreBenchmark

PATH = os.path.dirname(os.path.realpath(__file__))

ds = PDBBindRefined(root=PATH+'/data/pdbbind', name='test')
print(ds[0])

ds = TMScoreBenchmark(root=PATH+'/data/tmscore', name='test', use_precomputed=False)
print(ds.tm_score['1adr_']['1aho_'])
