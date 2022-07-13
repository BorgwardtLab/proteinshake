Quickstart
============

To get started, follow the `installation instructions <https://torch-pdb.readthedocs.io/en/latest/modules/datasets.html>`_.

Loading a dataset is one in one line:

.. code-block:: python

  from torch_pdb.datasets import GODataset
  proteins = GODataset(root="/tmp/var")

This will fetch the protein coordinates and annotations.

To use this dataset in a PyTorch model, we need to choose a way of encoding the raw coordinates into a data structure.

Here, we choose a residue graph, see `this page <https://torch-pdb.readthedocs.io/en/latest/modules/representations.html>`_ for a complete list of available data representations.

.. code-block:: python

   from torch_pdb.representations import GraphDataset

   graph_data = GraphDataset(root='/tmp/var', proteins=proteins.proteins)


