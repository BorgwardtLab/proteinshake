import numpy as np
from sklearn.model_selection import train_test_split

from proteinshake.datasets import ProteinLigandDecoysDataset
from proteinshake.tasks import Task

class VirtualScreenTask(Task):
    """ Test an affinity scoring model on a virtual screen.
    The goal in a virtual screen is: for a given protein and a library of potential
    binders, bring the binders to the top of the list.
    In this task, the model is given a protein and a list of ligands to score.
    The model scores each ligand in a library with a score proportional to the likelihood
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

    .. code-block:: python

        >>> from proteinshake.tasks import VirtualScreenTask
        >>> import numpy as np
        >>> task = VirtualScreenTask()
        # predict a (random) binding score for each molecule 
        >>> preds = [np.random.rand(len(task.target(p))) for p in task.dataset.proteins()]
        >>> task.evaluate(preds, cutoff_fraction=.2)
        {'enrichment_factor-@.2': 0.6}

    """

    DatasetClass = ProteinLigandDecoysDataset

    def __init__(self, *args, **kwargs):
        kwargs['split'] = 'none'
        super().__init__(*args, **kwargs)

    @property
    def task_type(self):
        return ('protein', 'virtual_screen')

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

    def evaluate(self, y_true, y_pred, cutoff_fraction=.2):
        """ Computing enrichment factor on the whole dataset.
        To compute the EF, for each active ligand we check if it is in the top
        ``cutoff_fraction`` of the scored ligands and if so, count it as a 1 and 0 otherwise.
        Taking the mean gives a score between 0 and 1 where 1 is a perfect score.

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
        
        efs = []
        for lig_ids, (i, protein) in zip(y_true, enumerate(self.proteins)):
            active_ids = lig_ids[:protein['protein']['num_ligands']]

            cutoff_index = int(len(lig_ids) * cutoff_fraction)

            scores_dict = {lig_ids[i]: score for i, score in enumerate(y_pred[i])}
            ranks_dict = {lig_id: i < cutoff_index for i, lig_id in  enumerate(sorted(lig_ids, key = lambda x: scores_dict[x]))}

            mean_active_rank = np.mean([ranks_dict[lig_id] for lig_id in active_ids])
            efs.append(mean_active_rank)

        return {'enrichment_factor-@{self.cutoff_fraction}': np.mean(efs)}

