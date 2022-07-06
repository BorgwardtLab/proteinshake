.. torch-pdb documentation master file, created by
   sphinx-quickstart on Fri Jul  1 16:33:46 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to torch-pdb!
=====================================

**torch-pdb** is a collection of **pytorch** datasets built from protein structure data banks.
Each dataset can be imported directly into your pyg models for training and inspection.
Depending on the dataset, specific annotations will be included in the `Data` objects such as small molecule binding pockets, and function annotations.

We currently draw from the following protein structure respositories:

* `RCSB <https://www.rcsb.org/>`_: central repository for 3D structures (contains Gene Ontology, Enzyme Classification, Protein Family annotations, etc.). 
* `TMAlign <https://zhanggroup.org/TM-align/>`_: curated benchmark dataset for protein-protein simlarity computations. 
* `AlphaFold <https://www.deepmind.com/open-source/alphafold-protein-structure-database>`_: database of predicted protein 3D structures. 

Drawing from these sources we provide the following datasets. Each dataset is defined by a set of proteins for which some property is known (e.g. Gene Ontology annotation):

+--------+-------+-----------+--------------------+-----+------------+
| Name   | N     | Mean Size | Property           | Val | Type       |
|        | umber | (#        |                    | ues |            |
|        | of    | residues) |                    |     |            |
|        | struc |           |                    |     |            |
|        | tures |           |                    |     |            |
+========+=======+===========+====================+=====+============+
| RCSBD  | 9117  | 421.945   | -                  | -   | -          |
| ataset |       |           |                    |     |            |
+--------+-------+-----------+--------------------+-----+------------+
| PfamD  | 7173  | 455.147   | Protein Family     | 2   | Ca         |
| ataset |       |           | (Pfam)             | 854 | tegorical, |
|        |       |           |                    | (ro | Hi         |
|        |       |           |                    | ot) | erarchical |
+--------+-------+-----------+--------------------+-----+------------+
| GOD    | 7811  | 442.024   | Gene Ontology (GO) | 73  | Ca         |
| ataset |       |           |                    | (ro | tegorical, |
|        |       |           |                    | ot) | Hi         |
|        |       |           |                    |     | erarchical |
+--------+-------+-----------+--------------------+-----+------------+
| ECD    | 1864  | 587.329   | Enzyme             | 633 | C          |
| ataset |       |           | Classification     |     | ategorical |
|        |       |           | (``EC``)           |     |            |
+--------+-------+-----------+--------------------+-----+------------+
| PD     | 5316  | 428.289   | Small Mol. Binding | 2   | Binary     |
| BBindR |       |           | Site               |     |            |
| efined |       |           | (residue-level)    |     |            |
+--------+-------+-----------+--------------------+-----+------------+
| TMSc   | 200   | 247.29    | TM Score           | [0  | Re         |
| oreBen |       |           |                    | -1] | al-valued, |
| chmark |       |           |                    |     | Pairwise   |
+--------+-------+-----------+--------------------+-----+------------+

.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Notes

   notes/installation


.. toctree::
   :glob:
   :maxdepth: 1
   :caption: Package Reference

   modules/datasets
   modules/utils


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
