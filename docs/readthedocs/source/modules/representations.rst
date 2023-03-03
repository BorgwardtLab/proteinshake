``proteinshake.representations``
===================================

Protein representations take a parsed protein structure file and convert it to a data structure which can be used by deep learning models.
We currently support graph, voxel, and point cloud.
These classes are called inside the dataset objects and you should not need to access them directly unless you want to add your own.


.. code-block:: python

   >>> from proteinshake.datasets import RCSBDataset
   >>> da = RCSBDataset().to_graph(eps=8)

.. currentmodule:: proteinshake.representations

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.representations.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.representations
    :members:
    :exclude-members: download, describe, add_protein_attributes, process, processed_file_names, raw_file_names, num_classes, get, get_raw_files, get_id_from_filename
