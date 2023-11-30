import unittest
import numpy as np
from proteinshake.frontend.dataset import Dataset
from proteinshake.frontend.task import Task
from proteinshake.frontend.transforms import *


class TestTask(unittest.TestCase):
    def test_task(self):
        # CREATOR
        class MyTargetTransform(TargetTransform):
            def __call__(self, protein):
                return protein["label"]

        class MyEvaluator:
            def __call__(self, y_true, y_pred):
                return {"Accuracy": np.random.random()}

        class MySplitter:
            def __call__(self, dataset):
                index = {"train": [], "test": [], "val": []}
                for i, p in enumerate(dataset.proteins):
                    index[p["split"]].append(i)
                return index

        class MyTask(Task):
            def __init__(self, **kwargs):
                super().__init__(
                    dataset=Dataset(),
                    splitter=MySplitter(),
                    target=MyTargetTransform(),
                    evaluator=MyEvaluator(),
                    **kwargs
                )

        # END USER
        y_transform = lambda label: label * 100
        X_transform = DataTransform(
            representation_transform=PointRepresentationTransform(),
            framework_transform=TorchFrameworkTransform(),
        )

        task = MyTask(
            X_transform=X_transform,
            y_transform=y_transform,
        )

        print(task.train_index)
        print(next(task.X_train).shape)
        for X, y in task.train_dataloader(batch_size=32):
            print("X", X.shape)
            print("y", y)
            break


if __name__ == "__main__":
    unittest.main()