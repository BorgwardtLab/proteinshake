
<p align="center">
<img src="https://github.com/BorgwardtLab/proteinshake/raw/main/docs/images/logo_subtitle.svg" width="70%">
</p>

[pypi-url]: https://pypi.org/project/proteinshake
![test workflow](https://github.com/BorgwardtLab/proteinshake/actions/workflows/build.yml/badge.svg)
[![Documentation Status](https://readthedocs.org/projects/proteinshake/badge/?version=latest)](https://proteinshake.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/proteinshake)](https://pypi.org/project/proteinshake/)
![PyPI - Downloads](https://img.shields.io/pypi/dm/proteinshake)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/LICENSE)
![visitors](https://visitor-badge.glitch.me/badge?page_id=BorgwardtLab.proteinshake)

* Fetch clean protein datasets in one line
* Convert proteins to graphs, point clouds, voxels, and surfaces (coming soon).
* Work in your favorite deep learning framework (pytorch, tensorflow, pytorch-geometric, dgl, networkx)



proteinshake is a collection of protein structure datasets built from [PDB](https://www.rcsb.org/) and [AlphaFold](https://alphafold.ebi.ac.uk/).
After installing, datasets can be passed directly to ML loaders for model training.

## Installation


```
$ pip install proteinshake
```

## Demo

How to load an AlphaFold dataset as of pytorch-geometric graphs.

```python
>>> from proteinshake.datasets import AlphaFoldDataset

>>> data = AlphaFoldDataset(root='.', organism='escherichia_coli').to_graph(k=5).pyg()
>>> protein_tensor, protein_data = data[0]
>>> protein_tensor
Data(x=[196], edge_index=[2, 0], edge_attr=[0, 1])
>>> protein_data['protein']['ID']
'P0A9H5'
>>> protein_data['protein']['sequence']
'MSDERYQQRQQRVKEKVDARVAQAQDERGIIIVFTGNGKGKTTAAFGTATRAVGHGKKVGVVQFIKGTWPNGERNLLEPHGVEFQVMATGFTWDTQNRESDTAACREVWQHAKRMLADSSLDMVLLDELTYMVAYDYLPLEEVVQALNERPHQQTVIITGRGCHRDILELADTVSELRPVKHAFDAGVKAQIGIDY'
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
