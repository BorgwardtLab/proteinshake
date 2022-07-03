![test workflow](https://github.com/BorgwardtLab/torch-pdb/actions/workflows/build.yml/badge.svg)

# `torch-pdb`: torch-geometric datasets built from the PDB

This is a collection of torch-geometric datasets built from [PDB](https://www.rcsb.org/).
After installing, datasets can be passed directly to torch loaders for model training.


## Installation


```
$ pip install torch-pdb
```

Note: ensure that you are using the correct versions of `torch-[scatter,sparse]` according to your hardware and cuda version. See [this](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html#installation-via-pip-wheels) page for more info.


### From source

```
$ git clone https://github.com/BorgwardtLab/torch-pdb
$ cd torch-pdb
$ pip install .
```

## Quickstart


```python
>>> from pdb_pyg.datasets import PDBBindRefined

>>> dataset = PDBBindRefined(name='pdbbind')
>>> d[0]
Data(edge_index=[2, 512], chain_id=[257], residue_idx=[257], residue=[257], residue_name=[257], residue_number=[257], residue_position=[257], coord=[257, 3], aa_idx=[257, 553], bond_type=[512], num_nodes=257, datapath='/tmp/var/test/raw/6ugp/6ugp_protein.pdb', name='6ugp')
```

## Datasets

This is a summary of available datasets.


| Dataset | # of graphs | graph attributes | node attributes | dataset attributes |
|---------|-------------|------------------|-----------------|-----------------|
| `RCSBDataset` | 9117 | | | |
| `GODataset` | 7811 | `GO` (`list`) | | |
| `ECDataset` | 1864 | `EC` (`list`) | | |
| `PfamDataset` | 7173 | `Pfam` (`list`) | | |
| `PDBBindRefined` | 5316 | | `is_site` (`binary`) | |
| `TMScoreBenchmark` | 200 | | | `tm_score[i][j]` (`float`), `rmsd[i][j]` (`float`) |

## Graph Building

You can customize the way graphs are built grom protein 3D structures with the following arguments to the dataset constructors:


*  `node_type (str)`: Currently only `residue` is supported.
*  `neighbor_type (str)`: We support `knn` and `radius`
*  `knn (int)`: Number of nearest neighbor residues to connect with an edge.


```python
from pdb_pyg.datasets import PDBBindRefined

dataset = PDBBindRefined(name='pdbbind', node_type='residue')
```

## Building dataset backends

If you want to compute some of the datasets we host (e.g. TM-score) follow these steps:

** [TM score benchmark dataset](https://zhanggroup.org/TM-align/benchmark/)

0. Install `TMalign` and add executable to `PATH` [download](https://zhanggroup.org/TM-align/)

1. Download the list of structures

```
$ python scripts/get_tmlist.py
```

2. Fetch PDBs and Compute TM scores

```
$  python scripts/tmscore_precompute.py --custom-urls torch_pdb/pkg_data/tm_pdblist.txt --dest data/tm-bench
```

3. Move the `tm-bench` folder to your `data` directory so that the `download()` method is skipped when creating the dataset.
