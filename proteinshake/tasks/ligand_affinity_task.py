class LigandAffinityTask(ShakeTask):
    def __init__(self, dataset, *args, **kwargs):
        super().__init__(dataset, *args, **kwargs)

    @property
    def target(self):
        return 'kd'

    def evaluator(self):
        return metrics.mean_squared_error

    def set_target_pyg(self, pyg_dataset, batch_size=32, framework='pyg'):
        d = []
        for p in pyg_dataset:
            p.y = getattr(p, self.target)
            d.append(p)
        return d


