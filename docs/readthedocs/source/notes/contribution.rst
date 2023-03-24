Contribution
============

We welcome any contributions (new tasks and datasets, features, or bug reports) through pull requests issues on the project `GitHub repository <https://github.com/BorgwardtLab/proteinshake>`_.

See the tutorial for :doc:`custom datasets and tasks<notebooks/custom.ipynb>` for an example on how to subclass ``Dataset`` and ``Task``.

Please write appropriate test cases for your code. New ``Dataset`` classes should be fully reproducible in our release pipeline, meaning that the raw PDB files need to be hosted on a publicly available server. New ``Task`` classes should implement appropriate metrics. The random, sequence, and structure splits will however be computed during the release, you don't have to take care of that.

If you have questions on how to contribute, please open an issue on GitHub or contact any of the authors.
