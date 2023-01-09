Protein Numerical Representations
===================================

.. currentmodule:: proteinshake.representations

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.representations.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.representations
    :members:
    :exclude-members: download, describe, add_protein_attributes, process, processed_file_names, raw_file_names, num_classes, get, get_raw_files, get_id_from_filename
