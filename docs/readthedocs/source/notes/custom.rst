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
        url = 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/UP000000805_243232_METJA_v4.tar'
        download_url(url, f'{self.root}/raw')
        extract_tar(f'{self.root}/raw/UP000000805_243232_METJA_v4.tar', f'{self.root}/raw/files')
        for path in glob.glob(f'{self.root}/raw/files/*.pdb.gz'):
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
        return glob.glob(f'{self.root}/raw/files/*.pdb')

ProteinShake also needs a unique ID to reference each individual protein, which we parse from the file name (in AlphaFoldDB files, this is the UniProt accession ID):

.. code:: python

    def get_id_from_filename(self, filename):
        return filename.rstrip('.pdb')

Lastly, the annotation step is implemented in the ``add_protein_attributes`` method.
Here we add the annotation to the ``protein_dict`` of each individual protein.

.. note::

    The ``protein_dict`` is the central storage item in ProteinShake.
    It contains the coordinates, meta data, and all annotations.
    See the ``Dataset`` source code documentation.

.. code:: python

    def add_protein_attributes(self, protein_dict):
        uniprot_id = self.protein_dict['protein']['ID']
        protein_dict['protein']['DNA-binding'] = self.annotations[uniprot_id]

.. tip::

    You can use the ``add_protein_attributes`` method also for filtering: if it returns ``None``, the protein will be removed from the dataset.

That's it! ProteinShake will now take care of downloading, parsing, cleaning and storing your data.
The whole code now looks like this:

.. code:: python

    import glob
    from proteinshake.dataset import Dataset
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
            return filename.rstrip('.pdb')

        def add_protein_attributes(self, protein_dict):
            uniprot_id = protein_dict['protein']['ID']
            dna_binding = self.annotations[uniprot_id]['ft_dna_bind']
            protein_dict['protein']['DNA-binding'] = dna_binding
            return protein_dict

Neat, right? You can use it like any other ProteinShake dataset:

.. code:: python

    dataset = DNABindingDataset()

.. tip::

    If you are happy with your custom dataset, consider :doc:`contributing it!<contribution>`


A custom task
-------------

