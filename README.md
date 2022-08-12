
<p align="center">
<img src="https://github.com/BorgwardtLab/proteinshake/raw/main/images/proteinshake_banner.png" width="70%">
</p>

# The largest repository of ML-ready protein 3D structure datasets

[pypi-url]: https://pypi.org/project/proteinshake
![test workflow](https://github.com/BorgwardtLab/proteinshake/actions/workflows/build.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/proteinshake/badge/?version=latest)](https://proteinshake.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/proteinshake)](https://pypi.org/project/proteinshake/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/proteinshake)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/LICENSE)
![visitors](https://visitor-badge.glitch.me/badge?page_id=BorgwardtLab.proteinshake)


This is a collection of protein structure datasets built from [PDB](https://www.rcsb.org/) and [AlphaFold](https://alphafold.ebi.ac.uk/).
After installing, datasets can be passed directly to ML loaders for model training.

## PDB Datasets
| name                                           |   num_proteins |   avg size (# residues) | property                                | values      | type                      |
|:-----------------------------------------------|---------------:|------------------------:|:----------------------------------------|:------------|:--------------------------|
| RCSBDataset                                    |          21989 |                 56.8898 | -                                     | -         | -                       |
| PfamDataset                                    |          18696 |                 59.4297 | Protein Family (Pfam)                   | 5215 (root) | Categorical, Hierarchical |
| GODataset                                      |          19267 |                 58.8485 | Gene Ontology (GO)                      | 101 (root)  | Categorical, Hierarchical |
| ECDataset                                      |           8150 |                 74.9618 | Enzyme Classification (`EC`)            | 2173        | Categorical               |
| PDBBindRefined                                 |           4642 |                108.806  | Small Mol. Binding Site (residue-level) | 2           | Binary                    |
| TMScoreBenchmark                               |            200 |                 49.458  | TM Score                                | [0-1]       | Real-valued, Pairwise     |

## AlphaFold Datasets
| name                                           |   num_proteins |   avg size (# residues) | property                                | values      | type                      |
|:-----------------------------------------------|---------------:|------------------------:|:----------------------------------------|:------------|:--------------------------|
| SwissProt          |          512.231 |                 79.334 | -                                     | -         | -                       |
| arabidopsis_thaliana          |          27434 |                 66.1312 | -                                     | -         | -                       |
| caenorhabditis_elegans        |          19694 |                 65.0678 | -                                     | -         | -                       |
| candida_albicans              |           5974 |                 62.782  | -                                     | -         | -                       |
| danio_rerio                   |          24664 |                 75.2797 | -                                     | -         | -                       |
| dictyostelium_discoideum      |          12622 |                 85.9275 | -                                     | -         | -                       |
| drosophila_melanogaster       |          13458 |                 81.2947 | -                                     | -         | -                       |
| escherichia_coli              |           4363 |                 51.5408 | -                                     | -         | -                       |
| glycine_max                   |          55799 |                 58.0664 | -                                     | -         | -                       |
| homo_sapiens                  |          23391 |                105.457  | -                                     | -         | -                       |
| methanocaldococcus_jannaschii |           1773 |                 46.7467 | -                                     | -         | -                       |
| mus_musculus                  |          21615 |                 83.0434 | -                                     | -         | -                       |
| oryza_sativa                  |          43649 |                 44.1931 | -                                     | -         | -                       |
| rattus_norvegicus             |          21272 |                 78.1547 | -                                     | -         | -                       |
| saccharomyces_cerevisiae      |           6040 |                 80.0745 | -                                     | -         | -                       |
| schizosaccharomyces_pombe     |           5128 |                 76.2427 | -                                     | -         | -                       |
| zea_mays                      |          39299 |                 46.1618 | -                                     | -         | -                       |



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

## Usage

See the [quickstart](https://proteinshake.readthedocs.io/en/latest/notes/quickstart.html) guide on our [documentation](https://proteinshake.readthedocs.io/en/latest/index.html) site to get started.

## Legal Note

<!---
We make our code available under the [MIT License](https://github.com/BorgwardtLab/proteinshake/blob/main/LICENSE). The datasets are distributed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).
-->

We obtained and modified data from the following sources:

The AlphaFold protein structures were downloaded from the [AlphaFold Structure Database](https://alphafold.ebi.ac.uk/), licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

The RCSB protein structures were downloaded from [RCSB](https://www.rcsb.org/), licensed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
