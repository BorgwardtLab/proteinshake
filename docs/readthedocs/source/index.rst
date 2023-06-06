.. torch-pdb documentation master file, created by
   sphinx-quickstart on Fri Jul  1 16:33:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ProteinShake!
=====================================

**ProteinShake** provides one-liner imports of large scale, preprocessed protein structure tasks and datasets for various model types and frameworks.

We provide a collection of preprocessed and cleaned protein 3D structure datasets from RCSB and AlphaFoldDB, including annotations. Structures are easily converted to graphs, voxels, or point clouds and loaded natively into PyTorch, TensorFlow, NumPy, JAX, PyTorch Geometric, DGL and NetworkX. The task API enables standardized benchmarking on a variety of tasks on protein and residue level.

Check out the :doc:`installation<notes/installation>`, the :doc:`quickstart guide<notes/quickstart>`, or `our website <https://borgwardtlab.github.io/proteinshake/>`_ to get started.

We welcome contributions and bug reports through issues and pull requests on `GitHub <https://github.com/BorgwardtLab/ProteinShake>`_.

Who is ProteinShake for?
-------------------------

ProteinShake is intended for computational biologists and machine learning resarchers who need quick access to ML-ready curated datasets. ProteinShake is a large repository of ready-to-use datasets with rich property annotations. Going beyond experimental structures, this is also making AlphaFold predicted structures readily available for large-scale deep learning models.

We put emphasis on extendability, aiming to eliminate boilerplate code for data preparation and model evaluation across machine learning disciplines. ProteinShake therefore also serves as a general framework for processing protein structure data, and we hope it will serve the community as a platform to share their datasets and evaluation tasks. New datasets and tasks can be created with just a few lines of code (see the :doc:`tutorial<notebooks/custom>`) and we will integrate your contributions through pull requests on GitHub (see the :doc:`contribution guide<notes/contribution>`) to be part of the ProteinShake releases. 

.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Notes

   notes/installation
   notes/quickstart
   notes/contribution
   notes/submission
   notes/citation

.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Tutorials

   notebooks/dataset
   notebooks/task
   notebooks/custom
     

.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Package Reference

   modules/datasets
   modules/tasks
   modules/representations
   modules/transforms
   modules/utils


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`