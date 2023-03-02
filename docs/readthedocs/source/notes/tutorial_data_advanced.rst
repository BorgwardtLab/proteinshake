Under the hood with datasets
-------------------------------

The core object is the ``Dataset``.
Constructing a dataset from scratch follows two major steps: downloading and parsing.
The result is an iterable collection of dictionaries which hold all the relevant information for the collection of proteins.
In the next section we provide more details on these processing steps.

Fetching Proteins
~~~~~~~~~~~~~~~~~~

Protein 3D data is conventionally stored in a coordinate file known as PDB, and more recently mmCIF.
In order to construct a dataset, proteinshake needs to obtain a list of PDB files for the proteins of interest.
Each dataset calls ``download()`` which places a PDB file for each protein of interest in `/path/to/root/raw/files/`.
In most cases this involves executing API requests to databases of protein files, but it can also simply copy over local PDB files if you have your own.

Most publicly available PDBs can be found in a handful of sources such as RCSB-PDB, AlphaFoldDB, and PDBBind so the downloading logic for those is already implemented by ProteinShake.
However you can still customize the downloading logic for these databses.
For example you can pass your own query to the RCSB download API: 

.. code-block:: python

        from proteinshake.datasets import RCSBDataset

        # fetches PDBs for which a Pfam annotation is known
        dataset = RCSBDataset(query=[['rcsb_polymer_entity_annotation.type','exact_match','Pfam'])

Or you can directly pass a list of PDBIDs of interest:

.. code-block:: python

        from proteinshake.datasets import RCSBDataset

        # fetches PDBs from a list of two IDs
        dataset = RCSBDataset(from_list=['5g04', '2hbs'])



Once constructed, every ``Dataset`` holds a list of protein structures which you can fetch by a simple list indexing:

.. code-block:: python

   >>> from proteinshake.Datasets import RCSBDataset
   >>> dataset = RCSBDataset(root="./data")
   # first protein
   >>> protein = dataset[0]


Default Protein Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each item in the `Dataset` is a multi-level dictionary which holds all the relevant information for the given protein.
Each dataset child class can add its own attributes, while transforms (LINK) can modify these attributes.
In the most basic datasets (``RCSBDataset`` and ``AlphaFoldDataset``) we store coordinates and chain information (if multi-chain), as well as surface accessibility values.

.. code-block:: python

        >>> protein['protein'] # protein-level data
        >>> protein['residue'] # residue-level data
        >>> protein['atom']    # atom-level data   

Each of the above three keys map to a dictionary which holds the actual annotations for each level of resolution.
Looking closer at the ``protein`` key we get the protein ID and its sequence:

.. code-block:: python

        >>> protein['protein']
        {'ID': '1JKW',
        'sequence': 'WTFSSEEQLARL...NRTCLSQLLDIMKSSEEVAVLKQKLDRCHSAELAL'}


The ``atom`` and ``residue`` keys map to dictionaries which hold lists of values the same size as the number of atoms and residues in the protein.

.. code-block:: python

        >>> protein['residue']
        {'residue_number': [11, 12, 13,...],
         'residue_type': ['W', 'T', 'F', 'S', 'S', 'E'..],
         'SASA': [...],
         'RSASA': [...],
         'x': [-19.697999954223633, -19.94300079345703, -18.11199951171875, ...]
         'y': [...],
         'z': [...],
        }

For atom resolution you have the same structure except the lists now have one entry per atom in the protein and instead of ``residue_type`` and ``residue_number`` we have ``atom_number`` and ``atom_type``.

Datasets for which all proteins have only a single chain(e.g. `RCBDataset` and `AlphaFoldDataset` there is no info on the chain). 
For others such as ``ProteinProteinInterfaceDataset`` where each protein is bound to another one, we have ``chain_id`` which tells you which chain each residue belongs to.


Custom Protein Attributes
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The two largest datasets (``RCSBDataset`` and ``AlphaFoldDataset``) only contain information about the protein structure itself.
However, we often collect information __about__ either the whole protein, or parts of the protein (residues, atoms, chains, substructures).
The ``Dataset`` object implements a method called ``add_protein_attributes(protein)`` which is applied to every raw protein dictionary at construction time and adds new keys to the dictionary.
The ``add_protein_attributes()`` can apply any logic to the protein such as computing surfaces, interfaces, or simply looking up annotations in a database of choice.

For example, the ``ProteinProteinInterfaceDataset`` has an additional key at the residue and atom levels which is ``True`` if the atom/residue is on the interface of two protein chains and ``False`` otherwise.
