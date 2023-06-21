
from typing import Any


class DummyModel():

    def __init__(self, task):
        self.task = task

    def train_step(self, *args, **kwargs):
        pass

    def test_step(self, *args, **kwargs):
        return self.task.dummy_output()

