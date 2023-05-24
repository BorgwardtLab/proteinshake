import numpy as np
from sklearn import metrics

from proteinshake.datasets import ProteinProteinInterfaceDataset
from proteinshake.tasks import Task
from proteinshake.transforms import CenterTransform, RandomRotateTransform, Compose

class ProteinProteinInterfaceTask(Task):
    """ Identify the binding residues of a protein-protein complex. This is a residue-level binary classification.

    .. code-block:: python

        >>> from proteinshake.tasks import ProteinProteinInterfaceTask
        >>> ta = ProteinProteinInterfaceTask()
    """

    DatasetClass = ProteinProteinInterfaceDataset

    @property
    def task_in(self):
        return ('residue', 'residue')

    @property
    def task_type(self):
        return ('residue', 'binary')

    @property
    def task_out(self):
        return ('binary')

    @property
    def out_dim(self):
        return (1)

    def dummy_output(self):
        import random
        return [random.randint(0, 1) for p in self.test_targets]

    def update_index(self):
        """ Transform to pairwise indexing """
        self.train_index = self.compute_pairs(self.train_index)
        self.val_index = self.compute_pairs(self.val_index)
        self.test_index = self.compute_pairs(self.test_index)
        
    def compute_targets(self):
        self.train_targets = [self.target(self.proteins[i], self.proteins[j]) for i,j in self.train_index]
        self.val_targets = [self.target(self.proteins[i], self.proteins[j]) for i,j in self.val_index]
        self.test_targets = [self.target(self.proteins[i], self.proteins[j]) for i,j in self.test_index]

    def compute_pairs(self, index):
        """ Grab all pairs of chains that share an interface"""
        def find_index(pdbid, chain):
            for i, p in enumerate(self.dataset.proteins()):
                if pdbid == p['protein']['ID'] and chain == p['residue']['chain_id'][0]:
                    return i
            raise IndexError

        proteins = self.dataset.proteins()
        chain_pairs = []
        for i, protein in enumerate(proteins):
            if i not in index:
                continue
            chain = protein['residue']['chain_id'][0]
            pdbid = protein['protein']['ID']
            try:
                chain_pairs.extend([(i, find_index(pdbid, partner)) for partner in self.dataset._interfaces[pdbid][chain]])
            # if chain is not in any interface, we skip
            except (KeyError, IndexError):
                continue
        return chain_pairs


    def target(self, protein_1, protein_2):
        chain_1 = protein_1['residue']['chain_id'][0]
        chain_2 = protein_2['residue']['chain_id'][0]
        chain_1_length = len(protein_1['residue']['chain_id'])
        chain_2_length = len(protein_2['residue']['chain_id'])
        pdbid = protein_1['protein']['ID']

        contacts = np.zeros((chain_1_length, chain_2_length))
        inds = np.array(self.dataset._interfaces[pdbid][chain_1][chain_2])
        contacts[inds[:,0], inds[:,1]] = 1.0
        return contacts

    @property
    def default_metric(self):
        return 'average_precision'

    def evaluate(self, y_true, y_pred):
        """ Evaluate performance of an interface classifier.
        """
        return {
            'auc_roc': metrics.roc_auc_score(y_true, y_pred),
            'average_precision': metrics.average_precision_score(y_true, y_pred),
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
