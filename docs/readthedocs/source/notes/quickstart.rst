Quickstart
==========

First, make sure you have :doc:`installed proteinshake<installation>`. 

Loading a dataset is a one-liner, consisting of three parts:

1. Choosing a dataset (e.g. AlphaFold or RCSB)

2. Converting proteins to a representation (graphs, point clouds, or voxels)

3. Loading them to your favorite framework (pyg, torch, dgl, tf, nx, np)

.. code-block:: python

  from proteinshake.datasets import RCSBDataset
  proteins = RCSBDataset(root="./data").to_point().torch()


The above line takes a dataset with proteins from the RCSB PDB Databank and transforms them to point clouds. They are then loaded into a PyTorch dataset (make sure you have PyTorch installed for this to work).

Some other examples:

.. code-block:: python

  # a graph dataset with epsilon-neighborhood graphs with radius 8 Angstrom, in DGL
  proteins = RCSBDataset(root="./data").to_graph(eps=8).dgl()

  # a graph dataset with k=5 nearest-neighbor graphs, in torch-geometric 
  proteins = RCSBDataset(root="./data").to_graph(k=5).pyg()

  # a voxel dataset with a voxel size of 10 Angstrom, in tensorflow
  proteins = RCSBDataset(root="./data").to_voxel(voxelsize=10).tf()

You can (almost) arbitrarily combine datasets, representations and frameworks.
These datasets can then be passed directly to the framework-specific dataloaders and models.


.. warning:: 

        Make sure you have installed the necessary downstream frameworks (e.g. torch, torch-geometric, etc.)
