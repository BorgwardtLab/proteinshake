How to contribute
=================

We welcome any contributions (new tasks and datasets, features, or bug reports) through pull requests on the project `GitHub repository <https://github.com/BorgwardtLab/proteinshake>`_.

After a review we will merge your pull request.
Your dataset/task will be processed on our server with the next release, and the pre-processed data will be available for download.
The random, sequence, and structure splits will be computed during the release, you don't have to take care of that.

If you have questions on how to contribute, please feel free to open an issue on GitHub or contact any of the authors.

.. tip::
    See the :doc:`Tutorial<custom>` for an example on how to build a custom ``Dataset`` and ``Task``.

.. important::

    The submitted dataset/task needs to fulfill the following requirements in order to be merged:

    - The dataset/task is well documented.
    - The dataset is fully reproducible and can be created anywhere (i.e. on our release server).
    - The dataset is compatible with our BSD-3/CC-BY-4.0 licenses (Please provide the source, license and, if applicable, citation in the docstring).
    - The task implements appropriate metrics.
    - Test cases have been implemented for the dataset/task.
    - The pull request passes the test pipeline (see the report on GitHub when opening a pull request. The tests will automatically run).

