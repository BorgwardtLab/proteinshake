Protein Datasets
========================

.. currentmodule:: proteinshake.datasets

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.datasets.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.datasets
    :members:
    :exclude-members: download, describe, add_protein_attributes, process, processed_file_names, raw_file_names, num_classes, get, get_raw_files, get_id_from_filename
