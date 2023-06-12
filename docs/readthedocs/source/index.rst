ProteinShake
=====================================

**ProteinShake provides one-liner imports of large scale, pre-processed protein structure datasets and tasks for various model types and frameworks.**

.. image:: /images/pipeline.png
  :width: 80%
  :align: center

We provide a collection of pre-processed and cleaned protein 3D structure datasets from RCSB and AlphaFoldDB, including annotations. Structures are easily converted to graphs, voxels, or point clouds and loaded natively into PyTorch, TensorFlow, NumPy, JAX, PyTorch Geometric, DGL and NetworkX. The task API enables standardized benchmarking on a variety of tasks on protein and residue level.

.. note::
   Check out the :doc:`Installation Guide<notes/installation>`, the :doc:`Quickstart<notes/quickstart>`, or our `Website <https://borgwardtlab.github.io/proteinshake/>`_ to get started.

.. note::
   We welcome contributions and bug reports through issues and pull requests on `GitHub <https://github.com/BorgwardtLab/ProteinShake>`_. See also our :doc:`Contribution Guide<notes/contribution>`.

Who is ProteinShake for?
-------------------------

ProteinShake is intended for computational biologists and machine learning researchers who need accessible datasets and well-defined evaluation benchmarks for their deep learning models.

We put emphasis on extendability, aiming to eliminate boilerplate code for data preparation and model evaluation across machine learning disciplines. ProteinShake therefore also serves as a general framework for processing protein structure data, and we hope it will serve the community as a platform to share their datasets and evaluation tasks. New datasets and tasks can be created with just a few lines of code (see the :doc:`Tutorial<notebooks/custom>`) and we will integrate your contributions through pull requests on GitHub (see the :doc:`Contribution Guide<notes/contribution>`). 

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
   modules/frameworks
   modules/transforms
   modules/utils