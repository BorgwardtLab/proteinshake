'''
Tests prediction task classes.
'''

import unittest, tempfile
from proteinshake.tasks import *

class TestTasks(unittest.TestCase):

    def task_check(self, task):
        self.assertIsNotNone(task.train_index)
        self.assertIsNotNone(task.train_targets)

    def test_ec_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = EnzymeClassTask(root=tmp, use_precomputed=False)
            self.task_check(task)

    def test_binding(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = BindingSiteDetectionTask(root=tmp, use_precomputed=False)
            self.task_check(task)

    def test_affinity(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = LigandAffinityTask(root=tmp, use_precomputed=False)
            self.task_check(task)

    def test_scop(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = StructuralClassTask(root=tmp, use_precomputed=False)
            self.task_check(task)

    def test_retrieve(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = StructureSimilarityTask(root=tmp, use_precomputed=False)
            self.task_check(task)


if __name__ == '__main__':
    unittest.main()
