import numpy as np
from sklearn import metrics

from proteinshake.datasets import ProteinProteinInterfaceDataset
from proteinshake.tasks import Task

class ProteinProteinInterfaceTask(Task):
    """ Identify the binding residues of a protein-protein complex.
    This is a residue-level binary classification.
    """

    DatasetClass = ProteinProteinInterfaceDataset
    
    type = 'Binary Classification'
    input = 'Protein and Protein'
    output = 'Protein Binding Interface Residues'
    default_metric = 'AUROC (median)'
    pairwise = True
    level = 'Residues'

    def compute_paired_index(self, index):
        """ Return all pairs of interactions that are in self.interfaces.
        """
        pdb_to_idx = {p['protein']['ID']:i for i,p in enumerate(self.dataset.proteins())}
        i,j = [],[]
        for k

    def target(self, protein_1, protein_2):
        pdbid_1, chain_1 = protein_1['protein']['ID'].split('_')
        pdbid_2, chain_2 = protein_2['protein']['ID'].split('_')
        shape = (len(protein_1['residue']['chain_id']), len(protein_2['residue']['chain_id']))
        return self.dataset.interfaces[pdbid_1][chain_1][chain_2], shape
    
    def target_transform(self, target):
        residue_indices, shape = target
        transformed_target = np.zeros(shape, dtype=bool)
        transformed_target[tuple(zip(*residue_indices))] = True
        return transformed_target

    def evaluate(self, y_true, y_pred):
        """ Evaluate performance of an interface classifier.
        """
        raw_values = {
            'auroc': np.zeros(len(y_true)),
            'auprc': np.zeros(len(y_true)),
            'sizes': np.zeros(len(y_true))
        }

        for i, (y, y_pred) in enumerate(zip(y_true, y_pred)):
            y = y.flatten()
            y_pred = y_pred.flatten()
            raw_values['auroc'][i] = metrics.roc_auc_score(y, y_pred)
            raw_values['auprc'][i] = metrics.average_precision_score(y, y_pred)
            raw_values['sizes'][i] = len(y)

        tot = np.sum(raw_values['sizes'])
        norm = raw_values['sizes'] / tot

        return {
            'AUROC (weighted)': np.mean(raw_values['auroc'] * norm),
            'AUPRC (weighted)':  np.mean(raw_values['auprc'] * norm),
            'AUROC (median)':  np.median(raw_values['auroc']),
            'AUPRC (median)':  np.median(raw_values['auprc']),
            'AUROC (mean)':  np.mean(raw_values['auroc']),
            'AUPRC (mean)':  np.mean(raw_values['auprc']),
        }
        
    def y_dummy(self):
        return [np.where(np.random.randint(0, 2, p.shape) == 0, 0, 1) for p in self.test_targets]