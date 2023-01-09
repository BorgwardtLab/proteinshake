
<p align="center">
<img src="docs/images/logo_subtitle.png#gh-light-mode-only" width="60%">
</p>
<p align="center">
<img src="docs/images/logo_subtitle_dark.png#gh-dark-mode-only" width="60%">
</p>

<div align="center">

![build](https://img.shields.io/github/actions/workflow/status/borgwardtlab/proteinshake/build.yml?color=%2303A9F4&style=for-the-badge)
[![pypi](https://img.shields.io/pypi/v/proteinshake?color=%2303A9F4&style=for-the-badge)](https://pypi.org/project/proteinshake/)
[![docs](https://img.shields.io/readthedocs/proteinshake?color=%2303A9F4&style=for-the-badge)](https://proteinshake.readthedocs.io/en/latest/?badge=latest)
[![downloads](https://img.shields.io/pypi/dm/proteinshake?color=%2303A9F4&style=for-the-badge)](https://pypi.org/project/proteinshake/)

</div>

<p align="center">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://borgwardtlab.github.io/proteinshake/#quickstart">Quickstart</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://borgwardtlab.github.io/proteinshake">Website</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://proteinshake.readthedocs.io/en/latest/?badge=latest">Documentation</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="">Paper</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://proteinshake.readthedocs.io/en/latest/notes/contributing.html">Contribute</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://borgwardtlab.github.io/proteinshake/#leaderboard">Leaderboard</a>
</p>

<div align="center">

### ProteinShake provides one-liner imports of large scale, preprocessed protein structure datasets for various model types and frameworks.

We provide a collection of preprocessed and cleaned protein 3D structure datasets from RCSB and AlphaFoldDB, including annotations. Structures are easily converted to graphs, voxels, or point clouds and loaded natively into PyTorch, Tensorflow, Numpy, JAX, PyTorch-Geometric, DGL and NetworkX. The task API enables standardized benchmarking on a variety of tasks on protein and residue level. Find more information on the <a href="https://borgwardtlab.github.io/proteinshake">Website</a> or the <a href="https://proteinshake.readthedocs.io/en/latest/?badge=latest">Documentation</a>.

**Installation:**

```
pip install proteinshake
```

</div>

---

**Example usage:**
```python
>>> from proteinshake.datasets import AlphaFoldDataset
>>> data = AlphaFoldDataset(organism='escherichia_coli').to_graph(k=5).pyg()
>>> graph, protein_dict = data[0]
>>> graph
Data(x=[196], edge_index=[2, 0], edge_attr=[0, 1])
>>> protein_dict['protein']['ID']
'P0A9H5'
>>> protein_dict['protein']['sequence']
'MSDERYQQRQQRVKEKVDARVAQAQDERGIIIVFTGNGK...'
```

<div align="center">

---

**Legal Note**

*We obtained and modified data from the following sources:*

*The AlphaFold protein structures were downloaded from the [AlphaFold Structure Database](https://alphafold.ebi.ac.uk/), licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).*

*The RCSB protein structures were downloaded from [RCSB](https://www.rcsb.org/), licensed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).*

</p>
