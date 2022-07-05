.. torch-pdb documentation master file, created by
   sphinx-quickstart on Fri Jul  1 16:33:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to torch-pdb!
=====================================

**torch-pdb** is a collection of **pytorch-geometric** datasets built from protein structure data banks.
Each dataset can be imported directly into your pyg models for training and inspection.
Depending on the dataset, specific annotations will be included in the `Data` objects such as small molecule binding pockets, and function annotations.


* PDBBind: contains binding site
* TMAlign: dataset used for benchmarking TMAlign protein-protein alignment tool 
* GO: each protein is labeled with its gene ontology  (GO) term
* RCSB-PDB: useful for unsupervised training
* AlphaFold: useful for unsupervised training



.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Package Reference

   modules/datasets
   modules/utils


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
