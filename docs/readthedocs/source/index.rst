.. torch-pdb documentation master file, created by
   sphinx-quickstart on Fri Jul  1 16:33:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ProteinShake!
=====================================

**ProteinShake** supports deep learning model development for protein 3D structures by providing:

* **Dataset loading** in one line
* **Multi-framework** support (PyTorch, Tensorflow, Numpy, JAX, PyTorch-Geometric, DGL, NetworkX)
* **Represent** proteins as graph, point cloud, and voxel.\
* **Large datasets** of experimental and computational protein structures
* **Prediction tasks** with easy splitting and evaluation
* **Customizable** tasks and datasets


Go to :doc:`quickstart<notes/quickstart>` guide to get started.

Source code is hosted on `GitHub <https://github.com/BorgwardtLab/ProteinShake>`_ and the project homepage is `here <https://borgwardtlab.github.io/proteinshake/>`_.

We welcome contributions and bug reports through issues and pull requests.

Who is ProteinShake for?
-------------------------

``ProteinShake`` is intended for computational biologists and machine learning resarchers who need quick access to ML-ready curated datasets. To our knowledge, this is the largest repository of stready to use datasets with rich property annotations. Going beyond experimental structures, this is also making AlphaFold predicted structures for SwissProt and over a dozen organism proteomes. 


.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Notes

   notes/installation
   notes/quickstart
   notes/datastats
   notes/contributing
   notes/citation

.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Tutorials

   notes/tutorial_data
   notes/tutorial_tasks
   notes/tutorial_data_advanced
     

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
