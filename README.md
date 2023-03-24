
<p align="center">
<img src="https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/docs/images/logo_subtitle.png#gh-light-mode-only" width="60%">
</p>
<p align="center">
<img src="https://raw.githubusercontent.com/BorgwardtLab/proteinshake/main/docs/images/logo_subtitle_dark.png#gh-dark-mode-only" width="60%">
</p>

<div align="center">

![build](https://img.shields.io/github/actions/workflow/status/borgwardtlab/proteinshake/build.yml?color=%2303A9F4&style=for-the-badge)
[![pypi](https://img.shields.io/pypi/v/proteinshake?color=%2303A9F4&style=for-the-badge)](https://pypi.org/project/proteinshake/)
[![docs](https://img.shields.io/readthedocs/proteinshake?color=%2303A9F4&style=for-the-badge)](https://proteinshake.readthedocs.io/en/latest/?badge=latest)
[![downloads](https://img.shields.io/pypi/dm/proteinshake?color=%2303A9F4&style=for-the-badge)](https://pypi.org/project/proteinshake/)
[![codecov](https://img.shields.io/codecov/c/gh/BorgwardtLab/proteinshake?color=%2303A9F4&style=for-the-badge&token=0NL6CQZ6MB)](https://codecov.io/gh/BorgwardtLab/proteinshake)
</div>

<p align="center">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://borgwardtlab.github.io/proteinshake/#quickstart">Quickstart</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://borgwardtlab.github.io/proteinshake">Website</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://proteinshake.readthedocs.io/en/latest/?badge=latest">Documentation</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="">Paper</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://proteinshake.readthedocs.io/en/latest/notes/contributing.html">Contribute</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://borgwardtlab.github.io/proteinshake/#leaderboard">Leaderboard</a>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
<a href="https://proteinshake.readthedocs.io/en/latest/notebooks/dataset.html">Tutorials</a>
</p>

<div align="center">

### ProteinShake provides one-liner imports of large scale, preprocessed protein structure tasks and datasets for various model types and frameworks.

We provide a collection of preprocessed and cleaned protein 3D structure datasets from RCSB and AlphaFoldDB, including annotations. Structures are easily converted to graphs, voxels, or point clouds and loaded natively into PyTorch, TensorFlow, NumPy, JAX, PyTorch Geometric, DGL and NetworkX. The task API enables standardized benchmarking on a variety of tasks on protein and residue level.

Find more information on the <a href="https://borgwardtlab.github.io/proteinshake">Website</a> and the <a href="https://proteinshake.readthedocs.io/en/latest/?badge=latest">Documentation</a>, or check out the <a href="https://proteinshake.readthedocs.io/en/latest/notebooks/dataset.html">Tutorials</a>.

**Installation:**
```diff
- This is a pre-release version. There may be unannounced changes to the API and datasets. -
- We expect some bugs as well, please open an issue if you find one. -
```
```
pip install proteinshake
```

</br>

**Data workflow:**

In one line you can import large datasets of protein 3D structures, encode them as graphs/voxel grids/point clouds, and port them to your favorite learning framework.

<p align="center">
  <img width="700" src="https://raw.githubusercontent.com/BorgwardtLab/proteinshake/anim/docs/images/data_animation.svg">
</p>
</p>

**Task workflow:**

The task API lets you easily access the underlying data for several tasks, get random/sequence/structure based splits, and evaluate your predictions.

<p align="center">
  <img width="700" src="https://raw.githubusercontent.com/BorgwardtLab/proteinshake/anim/docs/images/tasks_animation.svg">
</p>
</p>

</div>

---

<div align="center">

**Legal Note**

*Code in this repository is licensed under [BSD-3](https://github.com/BorgwardtLab/proteinshake/blob/main/LICENSE), the dataset files on Zenodo are licensed under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/legalcode).*

*We obtained and modified data from the following sources:*

*The AlphaFold protein structures were downloaded from the [AlphaFold Structure Database](https://alphafold.ebi.ac.uk/), available under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).*

*The RCSB protein structures were downloaded from [RCSB](https://www.rcsb.org/), available under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).*

*Protein and Ligand binding structures and annotations were downloaded from [PDBbind-CN](http://www.pdbbind.org.cn/), available under the [End User Agreement for Access to the PDBbind-CN Database and Web Site](http://www.pdbbind.org.cn/enroll.php).*

*The Gene Ontology was downloaded from the [Gene Ontology Consortium](http://geneontology.org/), available under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).*

*The SCOP data was downloaded from the [Structural Classification of Proteins](https://scop.mrc-lmb.cam.ac.uk/), available under [CC-BY-4.0](https://creativecommons.org/licenses/by/4.0/).*

</div>
