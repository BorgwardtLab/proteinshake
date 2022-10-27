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

    def __call(self, protein):
        protein.y = task.target(protein)
        return protein
