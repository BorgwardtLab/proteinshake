import numpy as np

from proteishake.datasets import ProteinLigandDecoysDataset
from proteinshake.tasks import Task

class VirtualScreenTask(Task):
    """ Test an affinity scoring model on a virtual screen

    .. warning::

        This is a zero-shot task so we use the whole dataset in
        evaluation. No train/test split.

    """
    DatasetClass = ProteinLigandDecoysDataset
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def compute_token_map(self):
        labels = {p['protein']['EC'].split(".")[self.ec_level] for p in self.proteins}
        return {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def task_type(self):
        return 'zero-shot'

    def target(self):
        """ The target here is a sorted list of smiles where the true ligands
        occupy the top of the list and the decoys occupy the rest.
        Since this is a zero-shot task we only use this internally for evaluation.

        """
        return protein['protein']['ligands_smiles'] + protein['protein']['decoys_smiles']

    @property
    def num_features(self):
        return 20

    def evaluate(self, y_pred, cutoff_fraction=.2):
        """ Computing enrichment factor on the whole dataset.
        To compute the EF, for each active ligand we check if it is in the top
        ``cutoff_fraction`` of the scored ligands and if so, count it as a 1 and 0 otherwise.
        Taking the mean gives a score between 0 and 1 where 1 is a perfect score.

        Arguments
        -----------

        y_pred: list[float]
            A list of binding scores for each molecule returned by ``self.target(protein)``.
            We assume that a large value for the score means stronger likelihood of binding.
        cutoff_fraction: float
            Top percentile based on given scores within which to count actives.


        Returns
        ----------
        metrics: dict
            A dictionary with the enrichment factor.

        """
        
        efs = []
        for protein in self.dataset.proteins():
            lig_ids = self.target(protein)
            active_ids = lig_ids[:protein['protein']['num_ligands']

            cutoff_index = int(len(lig_ids) * self.cutoff_percent)

            scores_dict = {lig_ids[i]: score for i, score in enumerate(y_pred)}
            ranks_dict = {lig_id: i < cutoff_index for i, lig_id in  sorted(lig_ids, key = lambda x: scores_dict[x])}

            mean_active_rank = np.mean([ranks_dict[lig_id] for lig_id in active_ids])
            efs.append(mean_active_rank)

        return {'enrichment_factor-@{self.cutoff_fraction}': np.mean(efs)}
