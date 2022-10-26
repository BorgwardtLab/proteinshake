from itertools import combinations

from sklearn import metrics
from sklearn.model_selection import train_test_split

from proteinshake.datasets import TMAlignDataset
from proteinshake.tasks import ShakeTask

class RetrieveTask(ShakeTask):
    def __init__(self, root, *args, **kwargs):
        dataset = TMAlignDataset(root=root)
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

    def target(self, idx):
        pdbid_1 = self.proteins[idx[0]]['protein']['ID']
        pdbid_2 = self.proteins[idx[1]]['protein']['ID']
        return self.dataset.tm_score[pdbid_1][pdbid_2]

    def evaluate(self, pred, true):
        return {'mse': metrics.mean_squared_error(pred, true)}

if __name__ == "__main__":
    # from proteinshake.datasets import TMAlignDataset
    # dataset = TMAlignDataset(root='tm')
    task = RetrieveTask('boo')
    print(task.target(task.train_ind[0]))
