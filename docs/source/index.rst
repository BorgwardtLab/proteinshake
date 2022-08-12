.. torch-pdb documentation master file, created by
   sphinx-quickstart on Fri Jul  1 16:33:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to proteinshake!
=====================================

With **proteinshake** you can quickly load datasets of 3D biomolecular structures into your favorite framework and structure representation. 

We currently collect structures from the following databases:

* `RCSB <https://www.rcsb.org/>`_: central repository for all solved 3D structures (contains Gene Ontology, Enzyme Classification, Protein Family annotations, etc.). 
* `TMAlign <https://zhanggroup.org/TM-align/>`_: curated benchmark dataset for protein-protein simlarity computations. 
* `PDBBind <http://www.pdbbind.org.cn/index.php/>`_: curated dataset of biomolecular complexes (protein-protein, protein-small molecule, etc.) including experimental affinity information. 
* `AlphaFold <https://www.deepmind.com/open-source/alphafold-protein-structure-database>`_: database of predicted protein 3D structures. 


Drawing from these sources we provide many large ML-ready `datasets <https://torch-pdb.readthedocs.io/en/latest/modules/datasets.html>`_ organized by annotation property or structure determination method.

Once a dataset is loaded, we provide several options for `representing <https://torch-pdb.readthedocs.io/en/latest/modules/representations.html>`_ the raw atomic coordinates into graphs, point clouds, and voxels (surfaces coming soon).

Go to `quickstart <https://torch-pdb.readthedocs.io/en/latest/notes/quickstart.html>`_ guide to get started.

Source code is hosted on `GitHub <https://github.com/BorgwardtLab/torch-pdb>`_.

We welcome contributions and bug reports through issues and pull requests.

Who is proteinshake for?
-------------------------

``proteinshake`` is intended for computational biologists and machine learning resarchers who need quick access to ML-ready curated datasets. To our knowledge, this is the largest repository of stready to use datasets with rich property annotations. Going beyond experimental structures, this is also making AlphaFold predicted structures for SwissProt and over a dozen organism proteomes. 


.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Notes

   notes/installation
   notes/quickstart
   notes/datastats
   notes/contributing


.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Package Reference

   modules/datasets
   modules/representations
   modules/utils


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
