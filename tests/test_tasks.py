'''
Tests prediction task classes.
'''

import random
import unittest, tempfile
from proteinshake.tasks import *


class TestTasks(unittest.TestCase):

    def task_check(self, task, skip_splits=False, pair=False):
        # check splits
        if not skip_splits:
            assert len(task.train_index) > 0
            assert len(task.test_index) > 0
            assert len(task.val_index) > 0
            assert len(set(task.train_index.flatten()).intersection(set(task.test_index.flatten()))) == 0
            assert len(set(task.train_index.flatten()).intersection(set(task.val_index.flatten()))) == 0
            assert len(set(task.val_index.flatten()).intersection(set(task.test_index.flatten()))) == 0

        # check targets
        prots = task.dataset.proteins()
        protein = next(prots)
        if pair:
            protein_2 = next(task.dataset.proteins())
            self.assertIsNotNone(task.target(protein, protein_2))
        else:
            self.assertIsNotNone(task.target(protein))

        # check eval
        self.assertIsNotNone(task.evaluate(task.test_targets, task.dummy_output()))

    def test_enzyme_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(EnzymeClassTask(split='random', root=tmp))
            self.task_check(EnzymeClassTask(split='sequence', root=tmp))
            self.task_check(EnzymeClassTask(split='structure', root=tmp))

    def test_protein_family(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(ProteinFamilyTask(split='random', root=tmp))
            self.task_check(ProteinFamilyTask(split='sequence', root=tmp))
            self.task_check(ProteinFamilyTask(split='structure', root=tmp))

    def test_gene_ontology(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(GeneOntologyTask(split='random', root=tmp))
            self.task_check(GeneOntologyTask(split='sequence', root=tmp))
            self.task_check(GeneOntologyTask(split='structure', root=tmp))

    def test_binding_site_detection(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(BindingSiteDetectionTask(split='random', root=tmp))
            self.task_check(BindingSiteDetectionTask(split='sequence', root=tmp))
            self.task_check(BindingSiteDetectionTask(split='structure', root=tmp))

    def test_ligand_affinity(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(LigandAffinityTask(split='random', root=tmp))
            self.task_check(LigandAffinityTask(split='sequence', root=tmp))
            self.task_check(LigandAffinityTask(split='structure', root=tmp))

    def test_protein_protein_interface(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(ProteinProteinInterfaceTask(split='random', root=tmp))
            self.task_check(ProteinProteinInterfaceTask(split='sequence', root=tmp))
            self.task_check(ProteinProteinInterfaceTask(split='structure', root=tmp))

    def test_structural_class(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(StructuralClassTask(split='random', root=tmp))
            self.task_check(StructuralClassTask(split='sequence', root=tmp))
            self.task_check(StructuralClassTask(split='structure', root=tmp))

    def test_structure_search(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(StructureSearchTask(split='random', root=tmp))
            self.task_check(StructureSearchTask(split='sequence', root=tmp))
            self.task_check(StructureSearchTask(split='structure', root=tmp))

    def test_virtual_screen(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(VirtualScreenTask(root=tmp), skip_splits=True)

    def test_structure_similarity(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.task_check(StructureSimilarityTask(split='random', root=tmp), pair=True)
            self.task_check(StructureSimilarityTask(split='sequence', root=tmp), pair=True)
            self.task_check(StructureSimilarityTask(split='structure', root=tmp), pair=True)

if __name__ == '__main__':
    unittest.main()
