
<div style="display:flex;justify-content:center;width:100%;">
<img src="docs/images/logo_subtitle.png" style="max-width:800px;width:100%;">
</div>

[![pypi](https://img.shields.io/pypi/v/proteinshake?color=%2303A9F4&style=for-the-badge)](https://pypi.org/project/proteinshake/)
![build](https://img.shields.io/github/actions/workflow/status/borgwardtlab/proteinshake/build.yml?color=%2303A9F4&style=for-the-badge)
[![docs](https://img.shields.io/readthedocs/proteinshake?color=%2303A9F4&style=for-the-badge)](https://proteinshake.readthedocs.io/en/latest/?badge=latest)
[![downloads](https://img.shields.io/pypi/dm/proteinshake?color=%2303A9F4&style=for-the-badge)](https://pypi.org/project/proteinshake/)
[![visitors]()]()

ProteinShake provides one-liner imports of large scale, preprocessed protein structure datasets for various model types and frameworks.

```
pip install proteinshake
```


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


<div style="display:flex;width:100%;align-items:center;justify-content:center;gap:20px;">
    <a href="" style="width:200px;height:40px;background-color:#03A9F4;color:#2d2d2d;text-align:center;">Quickstart</a>
    <a href="" style="width:200px;height:40px;background-color:#03A9F4;color:#2d2d2d;text-align:center;">Documentation</a>
    <a href="" style="width:200px;height:40px;background-color:#03A9F4;color:#2d2d2d;text-align:center;">Paper</a>
    <a href="" style="width:200px;height:40px;background-color:#03A9F4;color:#2d2d2d;text-align:center;">Contribute</a>
    <a href="" style="width:200px;height:40px;background-color:#03A9F4;color:#2d2d2d;text-align:center;">Leaderboard</a>
</div>

## Legal Note

<!---
We make our code available under the [BSD-3 License](https://github.com/BorgwardtLab/proteinshake/blob/main/LICENSE). The datasets are distributed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).
-->

We obtained and modified data from the following sources:

The AlphaFold protein structures were downloaded from the [AlphaFold Structure Database](https://alphafold.ebi.ac.uk/), licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).

The RCSB protein structures were downloaded from [RCSB](https://www.rcsb.org/), licensed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
