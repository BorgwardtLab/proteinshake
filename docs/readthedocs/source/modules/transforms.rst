``proteinshake.transforms``
=============================

Transforms are applied during the representation stage and you can pass them as a list to the conversion function. For example, this code snippet creates a dataset from RCSB and centers the coordinates of all the proteins::

        >>> from proteinshake.datasets import RCSBDataset
        >>> from proteinshake.transforms import CenterTransform
        >>> da = RCSBDataset().to_voxel(transforms=[CenterTransform()])


To create your own transform simply inherit from :meth:`proteinshake.transforms.Transform` and implement the ``__call__(self, protein)`` method.
This function is applied to the protein dictionary of every protein in a dataset.


.. currentmodule:: proteinshake.transforms

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.transforms.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.transforms
    :members:
    :exclude-members: target, evaluate, compute_splits, compute_token_map
