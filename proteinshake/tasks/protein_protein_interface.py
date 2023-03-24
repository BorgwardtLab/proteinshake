import numpy as np
from sklearn import metrics

from proteinshake.datasets import ProteinProteinInterfaceDataset
from proteinshake.tasks import Task
from proteinshake.transforms import CenterTransform, RandomRotateTransform, Compose

class ProteinProteinInterfaceTask(Task):
    """ Identify the binding residues of a protein-protein complex. This is a residue-level binary classification.

    NOTE: To make this an interesting task, the loader has to
    split the protein into its chains so that the model only sees
    one chain at a time. You should also pass the following transforms to the dataset
    :meth:`proteinshake.transforms.CenterTransform` and :meth:`proteinshake.transforms.RandomRotateTransform` when using representations that are aware of the atomic coordinates.

    NOTE: This task is currently in beta.


    .. code-block:: python

        >>> from proteinshake.tasks import ProteinProteinInterfaceTask
        >>> from proteinshake.transforms import CenterTransform, RandomRotateTransform
        >>> ta = ProteinProteinInterfaceTask()
        >>> data = ta.dataset.to_voxel(transforms=[CenterTransform(), RandomRotateTransform()).torch()
    """

    DatasetClass = ProteinProteinInterfaceDataset

    @property
    def task_type(self):
        return ('residue', 'binary') 

    def dummy_output(self):
        import random
        return [random.randint(0, 1) for p in self.test_targets]

    def compute_targets(self):
        # compute targets (e.g. for scaling)
        self.train_targets = [p for i in self.train_index for p in self.target(self.proteins[i])]
        self.val_targets = [p for i in self.val_index for p in self.target(self.proteins[i])]
        self.test_targets = [p for i in self.test_index for p in self.target(self.proteins[i])]

    def target(self, protein):
        return protein['residue']['is_interface']

    def evaluate(self, y_true, y_pred):
        """ Evaluate performance of an interface classifier.
        """
        return {
            'auc-roc': metrics.roc_auc_score(y_true, y_pred),
            'average precision': metrics.average_precision_score(y_true, y_pred),
        }

    def to_graph(self, *args, **kwargs):
        self.dataset = self.dataset.to_graph(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self

    def to_point(self, *args, **kwargs):
        self.dataset = self.dataset.to_point(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self

    def to_voxel(self, *args, **kwargs):
        self.dataset = self.dataset.to_voxel(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self
