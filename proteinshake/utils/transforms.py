class SetPyGTarget(object):
    """ Apply this to PyG dataset before model training.
    """
    def __init__(self, task):
        """ Initialize transform.

        Parameters
        -----------
        task: proteinshake.ShakeTask
            Proteinshake task object
        """

        self.task = task

    def __call__(self, data, protein_dict):
        data.y = self.task.target(protein_dict)
        return data, protein_dict
