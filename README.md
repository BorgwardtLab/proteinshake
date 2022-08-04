[pypi-url]: https://pypi.org/project/torch-pdb

![test workflow](https://github.com/BorgwardtLab/torch-pdb/actions/workflows/build.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/torch-pdb/badge/?version=latest)](https://torch-pdb.readthedocs.io/en/latest/?badge=latest)
[![PyPI version][pypi-img]][pypi-url]
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://raw.githubusercontent.com/BorgwardtLab/torch-pdb/main/LICENSE)
![visitors](https://visitor-badge.glitch.me/badge?page_id=BorgwardtLab.torch-pdb)

# `proteinshake`: the largest repository of ML-ready protein 3D structure datasets

This is a collection of torch-geometric datasets built from [PDB](https://www.rcsb.org/) and [AlphaFold](https://alphafold.ebi.ac.uk/).
After installing, datasets can be passed directly to ML loaders for model training.


| name                                           |   num_proteins |   avg size (# residues) | property                                | values      | type                      |
|:-----------------------------------------------|---------------:|------------------------:|:----------------------------------------|:------------|:--------------------------|
| RCSBDataset                                    |          21989 |                 56.8898 | -                                     | -         | -                       |
| PfamDataset                                    |          18696 |                 59.4297 | Protein Family (Pfam)                   | 5215 (root) | Categorical, Hierarchical |
| GODataset                                      |          19267 |                 58.8485 | Gene Ontology (GO)                      | 101 (root)  | Categorical, Hierarchical |
| ECDataset                                      |           8150 |                 74.9618 | Enzyme Classification (`EC`)            | 2173        | Categorical               |
| PDBBindRefined                                 |           4642 |                108.806  | Small Mol. Binding Site (residue-level) | 2           | Binary                    |
| TMScoreBenchmark                               |            200 |                 49.458  | TM Score                                | [0-1]       | Real-valued, Pairwise     |
| AlphaFoldDataset_arabidopsis_thaliana          |          27434 |                 66.1312 | -                                     | -         | -                       |
| AlphaFoldDataset_caenorhabditis_elegans        |          19694 |                 65.0678 | -                                     | -         | -                       |
| AlphaFoldDataset_candida_albicans              |           5974 |                 62.782  | -                                     | -         | -                       |
| AlphaFoldDataset_danio_rerio                   |          24664 |                 75.2797 | -                                     | -         | -                       |
| AlphaFoldDataset_dictyostelium_discoideum      |          12622 |                 85.9275 | -                                     | -         | -                       |
| AlphaFoldDataset_drosophila_melanogaster       |          13458 |                 81.2947 | -                                     | -         | -                       |
| AlphaFoldDataset_escherichia_coli              |           4363 |                 51.5408 | -                                     | -         | -                       |
| AlphaFoldDataset_glycine_max                   |          55799 |                 58.0664 | -                                     | -         | -                       |
| AlphaFoldDataset_homo_sapiens                  |          23391 |                105.457  | -                                     | -         | -                       |
| AlphaFoldDataset_methanocaldococcus_jannaschii |           1773 |                 46.7467 | -                                     | -         | -                       |
| AlphaFoldDataset_mus_musculus                  |          21615 |                 83.0434 | -                                     | -         | -                       |
| AlphaFoldDataset_oryza_sativa                  |          43649 |                 44.1931 | -                                     | -         | -                       |
| AlphaFoldDataset_rattus_norvegicus             |          21272 |                 78.1547 | -                                     | -         | -                       |
| AlphaFoldDataset_saccharomyces_cerevisiae      |           6040 |                 80.0745 | -                                     | -         | -                       |
| AlphaFoldDataset_schizosaccharomyces_pombe     |           5128 |                 76.2427 | -                                     | -         | -                       |
| AlphaFoldDataset_zea_mays                      |          39299 |                 46.1618 | -                                     | -         | -                       |



## Installation


```
$ pip install proteinshake
```

Note: ensure that you are using the correct versions of `torch-[scatter,sparse]` according to your hardware and cuda version. See [this](https://pytorch-geometric.readthedocs.io/en/latest/notes/installation.html#installation-via-pip-wheels) page for more info.


### From source

```
$ git clone https://github.com/BorgwardtLab/proteinshake
$ cd proteinshake
$ pip install .
```

### Usage

See the [quickstart](https://torch-pdb.readthedocs.io/en/latest/notes/quickstart.html) guide on our [documentation](https://torch-pdb.readthedocs.io/en/latest/index.html) site to get started.

### Licenses

We make our code available under the [MIT License](https://github.com/BorgwardtLab/torch-pdb/blob/main/LICENSE). The datasets are distributed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

We obtained and modified data from the following sources:

The AlphaFold protein structures were downloaded from the [AlphaFold Structure Database](https://alphafold.ebi.ac.uk/), licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

The RCSB protein structures were downloaded from [RCSB](https://www.rcsb.org/), licensed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
