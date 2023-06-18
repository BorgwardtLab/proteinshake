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

    def compute_pairs(self, index):
        # todo: would be easier to take the all-pairwise from the parent class and just filter if they are not in the interfaces
        """ Grab all pairs of chains that share an interface"""
        protein_to_index = {p['protein']['ID']: i for i, p in enumerate(self.dataset.proteins())}
        chain_pairs = []
        for i, protein in enumerate(self.dataset.proteins()):
            if i not in index: continue
            pdbid, chain = protein['protein']['ID'].split('_')
            try:
                chain_pairs.extend([
                    (i, protein_to_index[f'{pdbid}_{partner}'])
                    for partner in self.dataset._interfaces[pdbid][chain]
                ])
            except (KeyError, IndexError): # if chain is not in any interface, we skip
                continue
        # filter partner if not in same index (train/test/val)
        chain_pairs = [(i,j) for i,j in chain_pairs if i in index and j in index] # @carlos please check
        return tuple(map(tuple, zip(*chain_pairs)))

    def target(self, protein_1, protein_2):
        chain_1 = protein_1['residue']['chain_id'][0]
        chain_2 = protein_2['residue']['chain_id'][0]
        chain_1_length = len(protein_1['residue']['chain_id'])
        chain_2_length = len(protein_2['residue']['chain_id'])
        pdbid = protein_1['protein']['ID'].split('_')[0]

        print(pdbid, chain_1, chain_2)

        contacts = np.zeros((chain_1_length, chain_2_length))
        try:
            print(self.dataset._interfaces[pdbid][chain_1][chain_2])
            print('-----------')
            inds = np.array(self.dataset._interfaces[pdbid][chain_1][chain_2])
            contacts[inds[:,0], inds[:,1]] = 1.0
        except KeyError: # raised if there are no interactions between query chains
            pass
        return np.array(contacts)

    def evaluate(self, y_true, y_pred):
        """ Evaluate performance of an interface classifier.
        """
        raw_values = {
            'auroc': np.zeros(len(y_true)),
            'auprc': np.zeros(len(y_true)),
            'sizes':np.zeros(len(y_true))
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
        
    def dummy(self):
        return [np.where(np.random.randint(0, 2, p.shape) == 0, 0, 1) for p in self.test_targets]