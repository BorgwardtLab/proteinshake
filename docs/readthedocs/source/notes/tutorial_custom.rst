Custom datasets and tasks
----------------------------------------

Custom Dataset
~~~~~~~~~~~~~~~~

The ``Dataset`` object logic takes care of two major steps: (i) loading the raw PDBs, and (ii) parsing/annotating each protein.
Once these two are taken care of, all the downstream work (repsrentations, frameworks) can be re-used.
This makes customizing the datasets very simple.
Here is an example of how you can create a custom dataset in a situation where you have your own annotations on a local file that looks like this::

        pdbid, annotation
        xyz,2.3
        abc,1.0 cde,3.4


Each row corresponds to a protein that is hosted in the RCSB Databank so we can subclass the ``RCSBDataset`` object and add our own annotations::


        import pandas as pd
        from proteinshake.datasets import RCSBDataset

        class MyDataset(RCSBDataset):
                def __init__(self, annotation_csv, *args, **kwargs):
                        self.annotations = pd.read_csv(annotation_csv)
                        self.ids = self.annotations['pdbid']
                        super().__init__(from_list=self.ids, *args ,** kwargs)

                def add_protein_attribute(self, protein):
                        """ Store annotation in downloaded protein object"""
                        protein ["protein"]["my_annotation"] = self.annotations[protein.ID]["annotation"]
                        return protein


For more detailed information on how we construct datasets see (LINK).

Custom Task
~~~~~~~~~~~~

The ``Task`` object defines four major components: (i) the base dataset, (ii) the prediction target, (iii) train/validation/test splits, and (iv) prediction evaluation.

This is a template for a fully customized task::

        from proteinshake.tasks import Task
        from sklearn.metrics import mean_squared_error

        """
        (i) This static attribute needs to correspond to a
            proteinshake dataset from our hosted datasets
            or one that you created yourself
        """

        DatasetClass = MyTaskDataset

        class MyTask(Tasks):

                def target(self, protein):
                        """ (ii) Accepts a protein dictionary and
                            returns the prediction target for that protein.
                            This can accept pairs of proteins or other objects
                            depending on the task.
                        """ 
                        return protein['protein']['my_attribute']

                                
                def compute_custom_splits(self):
                        """ (iii) Return three lists of indices over 
                                  the task's dataset (``self.dataset``)
                        """
                        ...

                        return train_index, val_index, test_index

                def evaluate(self, y_pred):
                        """ (iv) Accepts a list of model outputs and returns 
                            a dictionary containing evaluation
                            metrics. By convention, ``y_pred``
                            is a list of values where each item corresponds
                            to a prediction on one item of the test set.
                        """
                        return {'mse': metrics.mean_squared_error(self.test_targets, y_pred)}

The user can now work with your task as follows::

        >>> import MyTask
        >>> task = MyTask()
        >>> data = task.dataset.to_graph(eps=8).pyg()
        >>> y_pred = model(data[task.train])
        >>> task.evaluate(y_pred)



.. tip ::
        
        If you would like to share your custom tasks and datasets, please feel free to open a `pull request <https://github.com/BorgwardtLab/proteinshake/pulls>`_ on our GitHub repository.



