class EnzymeCommissionTask(ShakeTask):
    def __init__(self, dataset, *args, **kwargs):

        super().__init__(dataset, *args, **kwargs)

    def compute_token_map(self):
        proteins, _ = self.dataset.proteins()
        labels = {self.label_process(p['protein'][self.target]) for p in proteins}
        self.token_map = {label: i for i, label in enumerate(sorted(list(labels)))}

    @property
    def target(self):
        return 'EC'

    def label_process(self, label):
        return  label.split(".")[0]

    @property
    def evaluator(self):
        return metrics.precision_score

    @property
    def level(self):
        return 'protein'


