from torch_pdb.datasets import PDBBindRefined
from torch_pdb.datasets import TMScoreBenchmark

d = PDBBindRefined(name='pdbind', root='data')
# d = TMScoreBenchmark(name='test', root='data')
# print(d.tm_scores[d[0].name][d[1].name])
