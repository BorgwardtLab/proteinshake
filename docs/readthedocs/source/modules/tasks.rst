``proteinshake.Task``
========================

These classes define the various prediction task we support.
See our :doc:`tasks tutorial <../notes/tutorial_tasks>` for a full example of the ``Task`` class, and :doc:`this tutorial <../notes/tutorial_custom>` to learn how to create your own tasks.

This snippet highlights the most common usage for a ``Task`` object. :: 

        >>> from proteinshake.tasks import EnzymeCommissionTask
        >>> dataset = task.dataset.to_graph(eps=8).pyg()
        >>> pred = model(dataset[task.train_index]) # assuming you have implemented model() elsewhere
        >>> task.evaluate(pred)
        {'precision': 0.5333515066547034, 'recall': 0.4799021029676011, 'accuracy': 0.6675514266755143}


.. currentmodule:: proteinshake.tasks

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.tasks.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.tasks
    :members:
    :exclude-members: target, evaluate, compute_splits, compute_token_map, task_type, num_classes, num_features, DatasetClass, compute_custom_split
