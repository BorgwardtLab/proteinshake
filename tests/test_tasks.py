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
            self.assertIsNotNone(task.train_index)
            self.assertIsNotNone(task.train_targets)

        # check targets
        prots = task.dataset.proteins()
        protein = next(prots)
        if pair:
            protein_2 = next(task.dataset.proteins())
            self.assertIsNotNone(task.target(protein, protein_2))
        else:
            self.assertIsNotNone(task.target(protein))

        # check eval
        pred = task.dummy_output()
        self.assertIsNotNone(task.evaluate(pred))

    def test_ec_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = EnzymeClassTask(split='structure', root=tmp)

    def test_pfam_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = ProteinFamilyTask(split='structure', root=tmp)

    # def test_go_task(self):
        # with tempfile.TemporaryDirectory() as tmp:
            # task = GeneOntologyTask(split='structure', root=tmp)
            # self.task_check(task)

    def test_binding(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = BindingSiteDetectionTask(split='structure', root=tmp)
            self.task_check(task)

    def test_affinity(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = LigandAffinityTask(split='structure', root=tmp)
            self.task_check(task)

    def test_ppi(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = ProteinProteinInterfaceTask(split='structure', root=tmp)
            self.task_check(task)

    def test_scop(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = StructuralClassTask(split='structure', root=tmp)
            self.task_check(task)

    def test_retrieve(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = StructureSearchTask(split='structure', root=tmp)
            self.task_check(task)

    def test_screen(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = VirtualScreenTask(split='structure', root=tmp)
            self.task_check(task, skip_splits=True)

    def test_sim(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = StructureSimilarityTask(split='structure', root=tmp)
            self.task_check(task, pair=True)

if __name__ == '__main__':
    unittest.main()
