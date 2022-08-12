Quickstart
============

To get started, follow the `installation instructions <https://torch-pdb.readthedocs.io/en/latest/modules/datasets.html>`_.

Loading a dataset is a one-liner, consisting of three parts:
1. Choosing a dataset (e.g. AlphaFold or RCSB)
2. Converting proteins to a representation (graphs, point clouds, or voxels)
3. Loading them to your favorite framework (pyg, torch, dgl, tf, nx, np)

.. code-block:: python

  from proteinshake.datasets import GODataset
  proteins = GODataset(root="./data").to_graph(k=5).pyg()

The above line takes a dataset with proteins annotated with Gene Ontology terms and transforms them to graphs with 5 neighbors. They are then loaded into a pytorch geometric dataset.

Some other examples:

.. code-block:: python

  # a graph dataset with epsilon-neighborhood graphs with radius 8 Angstrom, in DGL
  proteins = GODataset(root="./data").to_graph(eps=8).dgl()

  # a point cloud dataset, in tensorflow
  proteins = GODataset(root="./data").to_point().tf()

  # a voxel dataset with a voxel size of 10 Angstrom, in torch
  proteins = GODataset(root="./data").to_voxel(size=10).torch()

You can arbitrarily combine datasets, representations and frameworks.

The resulting dataset you can use like any other dataset of your framework. In the case of proteinshake for example:

.. code-block:: python

  from proteinshake.loader import DataLoader
  loader = DataLoader(proteins, batch_size=32)
