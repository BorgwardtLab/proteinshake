Quickstart
==========

.. note:: 

  If you haven't yet, :doc:`install proteinshake<installation>`. 

.. admonition:: **TLDR** (for the really impatient)

  .. code-block:: python

    dataset = EnzymeCommissionDataset()
    proteins = dataset.proteins(resolution='residue')
    protein_dict = next(proteins)

  .. code-block:: python

    task = EnzymeClassTask().to_point().torch()
    loader = DataLoader(task.train)
    metrics = task.evaluate(task.test_targets, my_predictions)

ProteinShake provides two main classes: the ``Dataset`` and the ``Task``.
Datasets are a collection of protein structures which can be converted to any of three representations (graphs, voxel grids, or point clouds) depending on your deep learning model.
Tasks extend the datasets with metrics and data splits for model evaluation.

Datasets
--------

Let's download a dataset from the ProteinShake data repository:

.. code-block:: python

  from proteinshake.datasets import EnzymeCommissionDataset
  dataset = EnzymeCommissionDataset()

This line will not do much however, because the files for the different resolutions (residue or atom) are lazily fetched from the data repository.
To start the download, we can manually access the dataset generator:

.. code-block:: python

  proteins = dataset.proteins(resolution='residue')

Now the download will start.
``dataset.proteins()`` returns a generator of dictionaries, each dictionary has two entries (``protein`` and ``residue``) which contain annotations on the protein and residue level, respectively.

We can for example access the sequence of the protein, which is a protein-level annotation:

.. code-block:: python

  protein_dict = next(proteins)
  sequence = protein_dict['protein']['sequence']

Other fields included in the ``protein_dict`` are the database ID, the coordinates (a residue-level property) and surface accessibility scores (RSA/SASA).
Each dataset also has fields specific to their proteins, such as the Enzyme Commission number (``protein_dict['protein']['EC']``) of the ``EnzymeCommissionDataset``.

Representations
---------------

The ``protein_dict`` is not perfectly suited for a deep learning model, which requires tensors of data.
Models may use different representations depending on their architecture.
Some models can work directly with the coordinates of the protein (the point cloud).
On the other hand, a 3D convolutional neural network could use a voxel grid representation, which rasterizes the protein on a regular grid.
And lastly, a graph neural network operates on a graph, which can be obtained by connecting residues with their neighbors.

ProteinShake therefore provides converters for all of these architectures.
They can simply be chained to the dataset which will result in a new dataset that is formatted accordingly:

.. code-block:: python
  :emphasize-text: .to_voxel()

  dataset = EnzymeCommissionDataset().to_voxel()

The ``.to_voxel()`` method can also be replaced with ``.to_graph()`` or ``.to_point()``.

There is not many things you can do with this dataset though, because it is a transitional object. Its purpose is to be converted further to your framework of choice:

.. caution:: 

  Make sure you have installed the necessary frameworks.

.. code-block:: python
  :emphasize-text: .torch()

  dataset = EnzymeCommissionDataset().to_voxel().torch()

where ``.torch()`` can be replaced by any of ``.tf()``, ``.np()``, ``dgl()``, ``.pyg()``, or ``.nx()``.

This object now is an iterable dataset which returns a tuple:

.. code-block:: python

  data, protein_dict = dataset[0]

The ``data`` object is the protein structure converted to (in this case) a voxel grid in PyTorch tensor format.
It depends on the converter and framework you used.
The ``protein_dict`` is the same as before.
You can use it to access annotations of the protein, for example to obtain class labels for your training.

To automate this, we provide a ``transforms`` interface.

.. code-block:: python

  def my_transform(item):
      data, protein_dict = item
      label = protein_dict['protein']['EC']
      # the EC number is a string looking like "1.3.5.14"
      # let's use only the first EC level as the label here
      label = int(label.split('.')[0])
      return data, torch.tensor(label)

  dataset = EnzymeCommissionDataset().to_voxel().torch(transform=my_transform)

This little snippet will automatically convert your target label to a tensor on the fly, everytime you access a protein.
The dataset ``__getitem__`` now returns a tuple with the protein structure data and a label tensor, which can conveniently be used during training.

The framework method also provides a ``pre_transform`` argument for applying the transform only once before saving the data, so it will not run every time you access a protein.
See the :doc:`Documentation<modules/frameworks>` for more information.


Tasks
-----

Now finally, of course we would like to use the datasets to do some model training and evaluation.
ProteinShake provides the ``Task`` classes, which extend the datasets with data splits and metrics.
They work very similar to a ``Dataset``:

.. code-block:: python

  from proteinshake.tasks import EnzymeClassTask
  task = EnzymeClassTask(split='sequence').to_voxel().torch()

.. note::

  You can change the ``split`` argument to retrieve either random, sequence, or structure-based splits.

The task has a few attributes and methods that are specific to model training and evaluation.
Let's look at our prediction targets.

.. code-block:: python

  print(task.test_targets)

We can retrieve the train, test and validation splits to put them into a dataloader.

.. note::

  ProteinShake is directly compatible with any dataloader from the supported frameworks. The usage may differ slightly. Check the `quickstart on the website <https://borgwardtlab.github.io/proteinshake/#quickstart>`_ for instructions.

.. code-block:: python

  from torch.utils.data import DataLoader
  train, test = DataLoader(task.train), DataLoader(task.test)

The task classes also implement appropriate metrics and function as an evaluator.

.. code-block:: python

  metrics = task.evaluate(task.test_targets, my_model_predictions)

This will return a dictionary of various relevant metrics.
Each task has a default metric which we use to rank models in the `Leaderboard <https://borgwardtlab.github.io/proteinshake/#leaderboard>`_.
Feel free to :doc:`submit your scores!<submission>`

.. code-block:: python

  print(metrics[task.default_metric])
