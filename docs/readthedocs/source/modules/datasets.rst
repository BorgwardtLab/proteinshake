``datasets``
==========================


These classes define our currently supported datasets. See our :doc:`datasets tutorial<../notes/dataset>` for a quick intro to datasets, and :doc:`this tutorial<../notes/custom>` to learn how to create your own datasets.

.. currentmodule:: proteinshake.datasets

.. autosummary::
    :nosignatures:
    {% for cls in proteinshake.datasets.classes %}
      {{ cls }}
    {% endfor %}

.. automodule:: proteinshake.datasets
    :members:
    :exclude-members: limit, describe, process, processed_file_names, raw_file_names, num_classes, get, get_raw_files, get_id_from_filename, check_arguments_same_as_hosted, download_precomputed, validate, download_complete, start_download
