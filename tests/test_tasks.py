'''
Tests prediction task classes.
'''

import unittest, tempfile
from proteinshake.tasks import (EnzymeCommissionTask,
                                BindingSitePredictionTask,
                                LigandAffinityTask,
                                SCOPTask,
                                RetrieveTask
                                )

class TestTasks(unittest.TestCase):

    def task_check(self, task):
        self.assertIsNotNone(task.train_ind)

    def test_ec_task(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = EnzymeCommissionTask(root=tmp)
            self.task_check(task)

    def test_binding(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = BindingSitePredictionTask(root=tmp)
            self.task_check(task)

    def test_affinity(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = LigandAffinityTask(root=tmp)
            self.task_check(task)

    def test_scop(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = SCOPTask(root=tmp)
            self.task_check(task)

    def test_retrieve(self):
        with tempfile.TemporaryDirectory() as tmp:
            task = RetrieveTask(root=tmp)
            self.task_check(task)


if __name__ == '__main__':
    unittest.main()
