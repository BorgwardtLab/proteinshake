``frameworks``
===================================

Datasets can natively be converted to various deep learning frameworks.
The base ``FrameworkDataset`` is an iterable dataset with a ``__getitem__`` method.
Each framework sub-class implements ``convert_to_framework`` and optionally ``load_transform`` to provide the framework specific conversion.


.. code-block:: python
  :emphasize-text: .pyg()

   from proteinshake.datasets import RCSBDataset
   dataset = RCSBDataset().to_graph(eps=8).pyg()

.. currentmodule:: proteinshake.frameworks

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.frameworks.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.frameworks
    :members:
    :exclude-members: convert_to_framework, load_transform
