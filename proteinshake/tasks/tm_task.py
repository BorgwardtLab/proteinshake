from itertools import combinations

from sklearn import metrics

class RetrieveTask(ShakeTask):
    def __init__(self, dataset, *args, **kwargs):

        super().__init__(dataset, *args, **kwargs)

    def compute_splits(self):
        print(f">>> computing splits")
        _, size = self.dataset.proteins()
        inds = list(combinations(range(size), 2))

        train, test = train_test_split(inds, test_size=1 - self.train_ratio)
        val, test = train_test_split(test, test_size=self.test_ratio/(self.test_ratio + self.validation_ratio))

        self.train_ind = train
        self.val_ind = val
        self.test_ind = test

    @property
    def target(self, idx):
        return self.dataset[idx[0]][idx[1]]['TM']

    @property
    def evaluator(self):
        return metrics.mean_squared_error
        pass

if __name__ == "__main__":
    from proteinshake.datasets import TMAlignDataset
    dataset = TMAlignDataset(root='tm')
    task = RetrieveTask(dataset)
