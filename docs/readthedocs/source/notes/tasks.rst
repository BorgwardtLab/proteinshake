Working with tasks
==================

.. note::

    Make sure you understand :doc:`how to work with datasets<datasets>` first.

Of course, at some point, we would like to use the datasets to train and evaluate a model.
ProteinShake provides the :class:`~proteinshake.tasks.Task` classes, which extend the datasets with data splits and metrics.
They work very similar to a :class:`~proteinshake.datasets.Dataset` in that they store a set of proteins with annotations, only with some additional functionality such as splits and evaluation methods:

.. code-block:: python

  from proteinshake.tasks import EnzymeClassTask
  task = EnzymeClassTask(split='sequence').to_voxel().torch()

.. hint::

  You can change the ``split`` argument to retrieve either ``'random'``, ``'sequence'``, or ``'structure'``-based splits.
  The latter two are based on sequence/structure similarity which we pre-compute for you.
  The split type influences how hard the generalization to the test set is for the model.

  The ``split_similarity_threshold`` argument controls the maximum similarity between train and test.

The task has a few attributes and methods that are specific to model training and evaluation.
Let's look at our prediction targets.

.. code-block:: python

  print(task.test_targets)

We can retrieve the train, test and validation splits to put them into a dataloader.

.. note::

  ProteinShake is directly compatible with any dataloader from the supported frameworks.
  The usage may differ slightly. Check the :doc:`Quickstart <quickstart>` to see the differences.

.. code-block:: python

  from torch.utils.data import DataLoader
  train, test = DataLoader(task.train), DataLoader(task.test)

The task classes also implement appropriate metrics and function as an evaluator.

.. tip::

  Every task implements a ``dummy_output()`` method you can use for testing if you don't have model predictions at hand. This method will return random values with the correct shape and type for the task.

  .. code:: python

    my_model_predictions = task.dummy_output()

.. code-block:: python

  metrics = task.evaluate(task.test_targets, my_model_predictions)

This will return a dictionary of various relevant metrics.
Each task has a default metric which we use to rank models in the `Leaderboard <https://borgwardtlab.github.io/proteinshake/#leaderboard>`_.
Feel free to :doc:`submit your scores!<submission>`

.. code-block:: python

  print(metrics[task.default_metric])
