import numpy as np
from sklearn.model_selection import train_test_split

from proteinshake.datasets import ProteinLigandDecoysDataset
from proteinshake.tasks import Task

class VirtualScreenTask(Task):
    """ Test an affinity scoring model on a virtual screen.
    The goal in a virtual screen is: for a given protein and a library of potential
    binders, bring the binders to the top of the list.
    In this task, the model is given a protein and a list of ligands to score.
    The model scores each ligand in a library with a score representing the likelihood
    that the protein and ligand will bind. This can be a docking score, energy calculation,
    or just a probability.
    Each protein's ligand library contains a certain number of active molecules (ligands)
    and a certai (larger) number of decoys (non-binders).
    We use the predicted scores to sort the whole library and calculate the position of each
    active ligand in the sorted library.
    Ligands in the topi percentiles which are known to be active contribute a 1 to the score
    and those below the cutoff contribute a 0.

    .. warning::
        This is a zero-shot task so we use the whole dataset in
        evaluation. No train/test split.

    .. admonition:: Task Card

        * **Input:** one protein, 
        * **Output:** list of molecules sorted by model
        * **Evaluation:** Enrichment Factor (Chen, Hongming, et al. "On evaluating molecular-docking methods for pose prediction and enrichment factors." Journal of chemical information and modeling 46.1 (2006): 401-415. Note: to keep scores between 0 and 1 we use a normalized version of EF which is the mean rank of active compounds  
)




    """

    DatasetClass = ProteinLigandDecoysDataset

    type = 'Ranking'
    input = 'Protein and Molecule'
    output = 'Affinity Score Ranking'

    def __init__(self, *args, **kwargs):
        kwargs['split'] = 'none'
        super().__init__(*args, **kwargs)
        self.test_targets = [self.target(p) for p in self.proteins]

    @property
    def task_in(self):
        return ('protein', 'molecule')

    @property
    def task_type(self):
        return ('protein', 'virtual_screen')

    @property
    def task_out(self):
        return ('virtual_screen')

    @property
    def target_dim(self):
        return (1)

    def target(self, protein):
        """ The target here is a sorted list of smiles where the true ligands
        occupy the top of the list and the decoys occupy the rest.
        Since this is a zero-shot task we only use this internally for evaluation.
        """
        return protein['protein']['ligands_smiles'] + protein['protein']['decoys_smiles']

    @property
    def num_features(self):
        return 20

    def dummy_output(self):
        import random
        return [[random.random() for _ in range(len(self.target(p)))] for p in self.proteins]

    @property
    def default_metric(self):
        return 'enrichment_factor'
        pass

    def evaluate(self, y_true, y_pred, cutoff_fraction=.2):
        """ Computing enrichment factor on the whole dataset.
        To compute the EF, we select a fraction of all ligands ``cutoff_fraction`` (typically 1-5 %), use the scoring function to sort the ligands and keep those in the top fraction.
        The EF is the ratio between two quantities: the fraction of the selection that contains active ligands, and the fraction of the whole dataset that contains actives. 
        Though not standard practice, we also report the mean-active rank which unlike EF is a value
        between 0 and 1 that reports the percentile of actives among the sorted ligand set.

        Arguments
        -----------
        y_pred: list[list[float]]
            A list of binding scores for each molecule returned by ``self.target(protein)``.
            We assume that a large value for the score means stronger likelihood of binding.
        cutoff_fraction: float
            Top fraction based on given scores within which to count actives.
        Returns
        ----------
        metrics: dict
            A dictionary with the enrichment factor.
        """

        efs,mars  = [], []
        for lig_ids, (i, protein) in zip(y_true, enumerate(self.proteins)):
            n_actives = protein['protein']['num_ligands']
            active_ids = lig_ids[:protein['protein']['num_ligands']]
            positions = list(range(len(active_ids)))

            is_active = [i < n_actives for i in range(len(lig_ids))] 
            screened = [status for _,status in sorted(zip(y_pred[i], is_active), reverse=True)]

            sel = int(len(lig_ids) * cutoff_fraction)
            n_actives_in_sel = sum(screened[:sel])
            mar = 1 - (np.mean([r for r, s in enumerate(screened) if s == 1]) / len(lig_ids))

            mars.append(mar)
            efs.append((n_actives_in_sel / sel ) / (n_actives / len(lig_ids)))

        return {'enrichment_factor': np.mean(efs),
                'mean_active_rank': np.mean(mars)}
