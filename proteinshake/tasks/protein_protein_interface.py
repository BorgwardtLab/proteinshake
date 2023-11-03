import numpy as np
from sklearn import metrics

from proteinshake.datasets import ProteinProteinInterfaceDataset
from proteinshake.tasks import Task
from proteinshake.transforms.transforms import Compose
from proteinshake.transforms.coords import CenterTransform, RandomRotateTransform

class ProteinProteinInterfaceTask(Task):
    """ Identify the binding interface of a protein-protein complex. Protein function is driven in large part by binding events between different protein chains to form 'complexes'. Understanding how proteins interact with each other has implications in unraveling complex biological mechanisms, and designing proteins with desirable interactions. The underlying data is taken from the PDBBind database. All pairs of residues belonging to different chains and coming from different protein chains within 6A of each other (Townshend et al., 2019)  are labeled as positive examples.

    .. admonition:: Task Summary 

        * **Input:** two protein chains
        * **Output:** binary label for each residue in both chains (1 if residue belongs to interface 0 otherwise)
        * **Evaluation:** AUROC (*Fout, Alex, et al. "Protein interface prediction using graph convolutional networks." Advances in neural information processing systems 30 (2017)*)


    """

    DatasetClass = ProteinProteinInterfaceDataset
    
    type = 'Binary Classification'
    input = 'Protein and Protein'
    output = 'Protein Binding Interface Residues'
    
    @property
    def num_classes(self):
        return 2

    @property
    def task_in(self):
        return ('residue', 'residue')

    @property
    def task_type(self):
        return ('residue_pair', 'binary')

    @property
    def task_out(self):
        return ('binary')

    @property
    def out_dim(self):
        return (1)

    def dummy_output(self):
        import random
        return [np.where(np.random.randint(0, 2, p.shape) == 0, 0, 1) for p in self.test_targets]

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
        protein_to_index = {p['protein']['ID']: i for i, p in enumerate(self.dataset.proteins())}
        def find_index(pdbid, chain):
            return protein_to_index[f'{pdbid}_{chain}']

        proteins = self.dataset.proteins()
        chain_pairs = []
        for i, protein in enumerate(proteins):
            if i not in index:
                continue
            #chain = protein['residue']['chain_id'][0]
            pdbid, chain = protein['protein']['ID'].split('_')
            try:
                chain_pairs.extend([(i, find_index(pdbid, partner)) for partner in self.dataset._interfaces[pdbid][chain]])
            # if chain is not in any interface, we skip
            except (KeyError, IndexError):
                continue
        chain_pairs = [(i,j) for i,j in chain_pairs if i in index and j in index] # @carlos please check
        return np.array(chain_pairs, dtype=int)


    def target(self, protein_1, protein_2):
        chain_1 = protein_1['residue']['chain_id'][0]
        chain_2 = protein_2['residue']['chain_id'][0]
        chain_1_length = len(protein_1['residue']['chain_id'])
        chain_2_length = len(protein_2['residue']['chain_id'])
        pdbid = protein_1['protein']['ID'].split('_')[0]

        contacts = np.zeros((chain_1_length, chain_2_length))
        try:
            inds = np.array(self.dataset._interfaces[pdbid][chain_1][chain_2])
            contacts[inds[:,0], inds[:,1]] = 1.0
        except KeyError: # raised if there are no interactions between query chains
            pass
        return np.array(contacts)

    @property
    def default_metric(self):
        return 'average_precision'

    def evaluate(self, y_true, y_pred):
        """ Evaluate performance of an interface classifier.
        """
        raw_values = {'auroc': np.zeros(len(y_true)),
                  'auprc': np.zeros(len(y_true)),
                  'sizes':np.zeros(len(y_true))
                   }

        for i, (y, y_pred) in enumerate(zip(y_true, y_pred)):
            y = y.flatten()
            y_pred = y_pred.flatten()
            raw_values['auroc'][i] = metrics.roc_auc_score(y, y_pred)
            raw_values['auprc'][i] = metrics.average_precision_score(y, y_pred)
            raw_values['sizes'][i] = len(y)

        result = {}
        tot = np.sum(raw_values['sizes'])
        norm = raw_values['sizes'] / tot
        result['auroc_weighted'] = np.mean(raw_values['auroc'] * norm)
        result['auprc_weighted'] = np.mean(raw_values['auprc'] * norm)
        result['auroc_median'] = np.median(raw_values['auroc'])
        result['auprc_median'] = np.median(raw_values['auprc'])
        result['auroc_mean'] = np.mean(raw_values['auroc'])
        result['auprc_mean'] = np.mean(raw_values['auprc'])

        return result

    def to_graph(self, *args, **kwargs):
        self.dataset = self.dataset.to_graph(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self

    def to_point(self, *args, **kwargs):
        self.dataset = self.dataset.to_point(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self

    def to_voxel(self, *args, **kwargs):
        self.dataset = self.dataset.to_voxel(*args, **kwargs, transform=Compose([CenterTransform(), RandomRotateTransform()]))
        return self
