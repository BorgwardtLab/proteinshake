Prediction Tasks
========================

.. currentmodule:: proteinshake.tasks

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.tasks.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.tasks
    :members:
    :exclude-members: target, evaluate, compute_splits, compute_token_map
