Custom datasets and tasks
=========================

ProteinShake is not only a repository for datasets and tasks, but also a framework for processing protein structure data.
Every dataset class implements all the steps necessary to reproduce the pre-processing of the data (before it gets hosted).
This makes it easy to derive new datasets and tasks, by re-using boilerplate code from ProteinShake.
When creating a new dataset or task, you only have to implement the processing logic specific to your application.

A custom dataset
----------------

For the sake of this tutorial we will implement a small dataset with AlphaFoldDB structures of *E. coli*.
We will query UniProtKB to obtain annotations for each protein and add them to the dataset.

Let's have a look at how the ProteinShake base ``Dataset`` class is structured.

.. code:: python

    class Dataset():

        def download(self):
            raise NotImplementedError

        def get_raw_files(self):
            raise NotImplementedError

        def get_id_from_filename(self, filename):
            raise NotImplementedError

        def add_protein_attributes(self, protein_dict):
            return protein_dict

.. note::

    We are omitting a few things here, have a look at the source code documentation if you are interested!


These are all the methods you need to implement when creating a custom dataset.
The most important one is the ``download`` method.
Here you implement the download of the raw protein data (and possibly annotation files).

For our example we download the data from AlphaFoldDB.
The aim is to have a collection of ``.pdb`` files in the dataset root directory.
By convention, we put these files into ``f'{self.root}/raw/files'``, although you can put them anywhere.

.. code:: python

    def download(self):
        base_url = 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/'
        file_name = 'UP000000625_83333_ECOLI_v4.tar'
        download_url(base_url+file_name, self.root+'/raw')
        extract_tar(self.root+'/raw/'+file_name, self.root+'/raw/files')
        for path in glob.glob(self.root+'/raw/files/*.pdb.gz'):
            unzip_file(path)

.. hint::

    ``download_url`` ``extract_tar`` ``unzip_file`` ``uniprot_query``

    We use a few convenience functions provided by ProteinShake here. They come from the ``utils`` module, check it out!

.. hint::

    ProteinShake provides general classes for AlphaFoldDB and RCSB PDB which implement the download logic.
    For this tutorial we write everything from scratch, but if your data comes from these databases we recommend basing your dataset on the provided classes.

Let's say we are interested in DNA binding and want to know whether a protein binds DNA or not.
To retrieve the annotations, we will query UniProtKB.

.. code:: python

    def download(self):
        ...
        self.annotations = uniprot_query('organism_id:83333', 'ft_dna_bind')

Here, ``organism_id:83333`` is a UniProt query that will return all *E.coli* proteins.
We request the ``ft_dna_bind`` column, which is the DNA-binding information.
The data is stored in a dictionary that is accessible via the UniProt ID.
We will use this later in the annotation step.

Next, we need to tell ProteinShake where to find the ``.pdb`` files.
We do this by implementing ``get_raw_files`` which returns a list of paths to each file.

.. code:: python

    def get_raw_files(self):
        return glob.glob(self.root+'/raw/files/*.pdb')

ProteinShake also needs a unique ID to reference each individual protein, which we parse from the file name (in AlphaFoldDB files, this is the UniProt accession ID):

.. code:: python

    def get_id_from_filename(self, filename):
        return filename.split('-')[1]

Lastly, the annotation step is implemented in the ``add_protein_attributes`` method.
Here we add the annotation to the ``protein_dict`` of each individual protein.

.. note::

    The ``protein_dict`` is the central storage item in ProteinShake.
    It contains the coordinates, meta data, and all annotations.
    See the ``Dataset`` source code documentation.

.. code:: python

    def add_protein_attributes(self, protein_dict):
        uniprot_id = protein_dict['protein']['ID']
        if not uniprot_id in self.annotations: return
        dna_binding = self.annotations[uniprot_id]['DNA binding']
        protein_dict['protein']['DNA binding'] = not dna_binding is None
        return protein_dict

.. tip::

    You can use the ``add_protein_attributes`` method for filtering: if it returns ``None``, the protein will be removed from the dataset.

That's it! ProteinShake will now take care of downloading, parsing, cleaning and storing your data.
The whole code now looks like this:

.. code:: python

    import glob
    from proteinshake.datasets import Dataset
    from proteinshake.utils import *

    class DNABindingDataset(Dataset):

        def download(self):
            base_url = 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/'
            file_name = 'UP000000625_83333_ECOLI_v4.tar'
            download_url(base_url+file_name, self.root+'/raw')
            extract_tar(self.root+'/raw/'+file_name, self.root+'/raw/files')
            for path in glob.glob(self.root+'/raw/files/*.pdb.gz'):
                unzip_file(path)
            self.annotations = uniprot_query('organism_id:83333', 'ft_dna_bind')

        def get_raw_files(self):
            return glob.glob(self.root+'/raw/files/*.pdb')

        def get_id_from_filename(self, filename):
            return filename.split('-')[1]

        def add_protein_attributes(self, protein_dict):
            uniprot_id = protein_dict['protein']['ID']
            if not uniprot_id in self.annotations: return
            dna_binding = self.annotations[uniprot_id]['DNA binding']
            protein_dict['protein']['DNA binding'] = not dna_binding is None
            return protein_dict

Neat, right? You can use it like any other ProteinShake dataset:

.. code:: python

    dataset = DNABindingDataset()

.. tip::

    If you are happy with your custom dataset, consider :doc:`contributing it!<contribution>`


A custom task
-------------

A dataset becomes truly valuable when you define how to evaluate a model on it.
In ProteinShake, this is called a task.
It comprises train/test/validation splits and metrics that assess the performance of the model.
The metrics depend on the label(s) that you are interested in.

We will create a task based on our custom ``DNABindingDataset``.
An empty task looks like this:

.. code:: python

    class Task:

        DatasetClass = None
        type = None
        input = None
        output = None

        def target(self, protein):
            raise NotImplementedError

        def evaluate(self, y_true, y_pred):
            raise NotImplementedError

First we need to tell ProteinShake which dataset this task is based on.
For this we assign the ``DatasetClass`` class attribute:

.. code:: python

    class DNABindingTask(Task):
        DatasetClass = DNABindingDataset

Then there are a few key properties that define how a task is structured.
The properties are ``type``, ``input`` and ``output``.
Models can query these attributes to make task-specific decisions, such as the number of output neurons, or the type of loss to be used.

.. code:: python

    class DNABindingTask(Task):
        ...
        type = 'Binary Classification'
        input = 'Protein'
        output = 'DNA Binding'

.. note::

    The ``type`` and ``input`` attribute have to follow a convention.
    See the task documentation for details.

The most important methods of a task are ``target`` and ``evaluate``.
The first defines how the prediction target value can be read from the ``protein_dict``, the latter defines a dictionary of appropriate metrics.
Let's implement the two.

.. code:: python

    def target(self, protein_dict):
        return protein_dict['protein']['DNA binding']

    def evaluate(self, y_true, y_pred):
        return {
            'Accuracy': sklearn.metrics.accuracy_score(y_true, y_pred),
            'MCC': sklearn.metrics.matthews_corrcoef(y_true, y_pred),
        }

.. tip::

    By default, a random split will be computed on the fly when you use the task.
    You can implement ``compute_custom_split`` to define your own splitting logic.

    The random, sequence, and structure splits will only be computed during a release.
    If you :doc:`contribute your task<contribution>` we will compute and host them for you.

And we are done with the task!
The whole class looks like the following.
Again, you can use it like any other ProteinShake task, convert them to a repesentation, and load them to your favorite framework dataloader.

.. code:: python

    import sklearn
    from proteinshake.tasks import Task

    class DNABindingTask(Task):

        DatasetClass = DNABindingDataset
        type = 'Binary Classification'
        input = 'Protein'
        output = 'DNA Binding'

        def target(self, protein_dict):
            return protein_dict['protein']['DNA binding']

        def evaluate(self, y_true, y_pred):
            return {
                'Accuracy': sklearn.metrics.accuracy_score(y_true, y_pred),
                'MCC': sklearn.metrics.matthews_corrcoef(y_true, y_pred),
            }
