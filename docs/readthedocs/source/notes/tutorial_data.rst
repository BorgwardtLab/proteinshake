Getting started with protein datasets
----------------------------------------

ProteinShake implements all the steps from raw protein coordinate files (PDB/mmCIF) to training-ready dataset.
We also host the result of these computations for several datasets which are listed :doc:`here<datasets>`.

You can obtain a dataset ready for model training in one line. 

Here is an example of how you would get a ``torch-geometric.Dataset`` object of proteins pulled from the RCSB PDB Data Bank as epsilon graphs.

.. code-block:: python

  from proteinshake.datasets import RCSBDataset

  # a graph dataset with epsilon-neighborhood graphs with radius 8 Angstrom, in DGL
  proteins = RCSBDataset(root="./data").to_graph(eps=8).pyg()

The above code executes the three main steps of dataset preparation:

1. Loading the processed protein data:  ``RCSBDataset(root="data")``
2. Converting the proteins to your representation of choice: ``.to_graph(eps=8)``
3. Converting the dataset to your framework of choice: ``.pyg()``

To reproduce the processing you can pass the ``use_precomputed=False`` flag to ``RCSBDataset()``.
This executes all the processing steps locally, whereas by default we try to fetch the dataset from the datasets we host on `Zenodo <https://sandbox.zenodo.org/record/1170307>`_  as the processing can be quite time-consuming.

Next, we break down the three steps into some more detail.

Loading protein data
~~~~~~~~~~~~~~~~~~~~~~~

The first step in the snippet above does most of the leg work.
Once the dataset object is created, it holds an iterable of dictionaries, one for each protein in the dataset which is accessed through the ``Dataset.proteins()`` method.

.. code-block:: python

        >>> from proteinshake.datasets import RCSBDataset
        >>> dataset = RCSBDataset(root="./data")
        >>> proteins = dataset.proteins()
        >>> next(proteins)
        {'protein': {...},
         'residue': {...}
        }

Different implementations of the ``Dataset`` parent class let you customize this step.
For example, the ``RCSBDataset`` accepts a ``from_list`` argument which lets you specify which PDBs to fetch.

.. warning:: 

   The default behaviour is to fetch the pre-processed proteins from Zenodo. In order to see custom parameters in the dataset creation stage, set ``use_precomputed=False``.



Protein representations
~~~~~~~~~~~~~~~~~~~~~~~~

Once the processed protein data is in the ``Dataset`` object we need to convert it to a representation that works with downstream deep learning models.
We currently support graphs, point clouds, and voxels.

.. code-block:: python

        >>> point_dataset = dataset.to_point()
        >>> graph_dataset = dataset.to_graph(eps=9)
        >>> voxel_dataset = dataset.to_voxel()


These methods can be applied to any ``proteinshake.Dataset`` child class and perform the necessary computations for the different representations.
Some notes on each representation:


Point Clouds
__________________

Point clouds simply return the x, y, z coordinates of the alpha carbon for each residue.
You can finds some relevant processing in the transforms module (LINK) such as centering and rotating.
For example, to center the point clouds:

.. code-block:: python

        >>> from proteinshake.transforms import CenterTransform
        >>> dataset = RCSBDataset(root="./data_transformed", transforms=[CenterTransform()])


Graphs
________

Graph construction can be done in two ways: epsilon or k nearest neighbors.
The epsilon graph is chosen when ``eps`` is passed with a distance threshold.
All pairs of residues within the distance threshold are connected by an edge.
If ``k=4`` is passed then each residue is connected by an edge to its 4 nearest neighbors.::

        >>> knn_graph = dataset.to_graph(k=4)
        >>> eps_graph = dataset.to_graph(eps=8)


You can obtain a weighted graph where weights correspond to the distance between connected residues::

        >>> eps_graph = dataset.to_graph(eps=8, weighted=True)


Voxels
________

For the voxel representation we place a 3D grid of voxels over the protein and include a one-hot encoding of the amino acid or atom types present at the each voxel. 

Frameworks
~~~~~~~~~~~~~~

The final step is converting the protein representation to a computation framework of choice (e.g. pytorch-geometric, dgl, JAX/Numpy, etc.)
frameworks are available for each task and that is how we end up with the complete dataset creation command: ::

        >>> data = dataset.to_graph(eps=8, weighted=True).pyg()


In this example we converted to pytorch-geometric objects but you can use many others. See the :doc:`Representation <../modules/representations>` page for more.
At this point you can pass the dataset to a dataloader in your framework of choice.



