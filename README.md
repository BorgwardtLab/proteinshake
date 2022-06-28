# `torch-pdb`: torch-geometric datasets built from the PDB

This is a collection of torch-geometric datasets built from [PDB](https://www.rcsb.org/). 
After installing, datasets can be passed directly to torch loaders for model training.


## Installation


```
$ pip install torch-pdb
```
		

### From source

```
git clone https://github.com/BorgwardtLab/torch-pdb
cd torch-pdb
pip install .
```

## Quickstart


```
from pdb_pyg.datasets import PDBBindRefined

dataset = PDBBindRefined(name='pdbbind')
```

## Datasets

This is a summary of available datasets.


| Dataset | # of graphs | graph attributes | node attributes | edge attributes |
|---------|-------------|------------------|-----------------|-----------------|
| `PDBBindRefined`         |  5316           |     `is_site (binary)`             |                 |                 |

## Graph Building

You can customize the way graphs are built grom protein 3D structures with the following arguments to the dataset constructors:


*  `node_type (str)`: Currently only `residue` is supported. 
*  `neighbor_type (str)`: We support `knn` and `radius`
*  `knn (int)`: Number of nearest neighbor residues to connect with an edge.


```
from pdb_pyg.datasets import PDBBindRefined

dataset = PDBBindRefined(name='pdbbind', node_type='residue')
```
