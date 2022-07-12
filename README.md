[pypi-img]: https://img.shields.io/pypi/v/torch-pdb

[pypi-url]: https://pypi.org/project/torch-pdb

![test workflow](https://github.com/BorgwardtLab/torch-pdb/actions/workflows/build.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/torch-pdb/badge/?version=latest)](https://torch-pdb.readthedocs.io/en/latest/?badge=latest)
[![PyPI version][pypi-img]][pypi-url]
![visitors](https://visitor-badge.glitch.me/badge?page_id=BorgwardtLab.torch-pdb)


# `torch-pdb`: torch-geometric datasets built from the PDB

![](images/torch-pdb.png)

This is a collection of torch-geometric datasets built from [PDB](https://www.rcsb.org/).
After installing, datasets can be passed directly to torch loaders for model training.


| name             |   num_proteins |   avg size (# residues) | property                                | values      | type                      |
|:-----------------|---------------:|------------------------:|:----------------------------------------|:------------|:--------------------------|
| RCSBDataset      |           9117 |                 421.945 | nan                                     | nan         | nan                       |
| PfamDataset      |           7173 |                 455.147 | Protein Family (Pfam)                   | 2854 (root) | Categorical, Hierarchical |
| GODataset        |           7811 |                 442.024 | Gene Ontology (GO)                      | 73 (root)   | Categorical, Hierarchical |
| ECDataset        |           1864 |                 587.329 | Enzyme Classification (`EC`)            | 633         | Categorical               |
| PDBBindRefined   |           5316 |                 428.289 | Small Mol. Binding Site (residue-level) | 2           | Binary                    |
| TMScoreBenchmark |            200 |                 247.29  | TM Score                                | [0-1]       | Real-valued, Pairwise     |



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


## Graph Building

You can customize the way graphs are built grom protein 3D structures with the following arguments to the dataset constructors:


*  `node_type (str)`: Currently only `residue` is supported.
*  `neighbor_type (str)`: We support `knn` and `radius`
*  `knn (int)`: Number of nearest neighbor residues to connect with an edge.


```python
from pdb_pyg.datasets import PDBBindRefined

dataset = PDBBindRefined(name='pdbbind', node_type='residue')
```
