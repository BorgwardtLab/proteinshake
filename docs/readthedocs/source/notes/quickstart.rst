Quickstart
============

First, make sure you have :doc:`installed proteinshake<installation>`. 

Loading a dataset is a one-liner, consisting of three parts:
1. Choosing a dataset (e.g. AlphaFold or RCSB)
2. Converting proteins to a representation (graphs, point clouds, or voxels)
3. Loading them to your favorite framework (pyg, torch, dgl, tf, nx, np)

.. code-block:: python

  from proteinshake.datasets import RCSBDataset
  proteins = RCSBDataset(root="./data").to_point(size=10).torch()


The above line takes a dataset with proteins from the RCSB PDB Databank and transforms them to graphs with 5 neighbors. They are then loaded into a pytorch geometric dataset (make sure you have pytorch-geometric installed for this to work).

Some other examples:

.. code-block:: python

  # a graph dataset with epsilon-neighborhood graphs with radius 8 Angstrom, in DGL
  proteins = RCSBDataset(root="./data").to_graph(eps=8).dgl()

  # a graph dataset, in torch-geometric 
  proteins = RCSBDataset(root="./data").to_graph().pyg()

  # a voxel dataset with a voxel size of 10 Angstrom, in torch
  proteins = RCSBDataset(root="./data").to_voxel(size=10).torch()

You can arbitrarily combine datasets, representations and frameworks.
These datasets can then be passed directly to framework-specific dataloaders and models.


.. warning:: 

        Make sure you have installed the necessary downstream frameworks (e.g. pytorch, torch-geometric, etc.)
