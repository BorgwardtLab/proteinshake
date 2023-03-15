import json

from proteinshake.datasets import RCSBDataset

class ProteinFamilyDataset(RCSBDataset):
    """ Proteins from RCSB for which the protein family (Pfam) term is known.
    Each protein in the dataset has a `Pfam` attribute which stores the list of protein families.

    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 31109


   .. list-table:: Annotations
      :widths: 25 35 45
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Pfam accession code
        - :code:`protein['protein']['Pfam']`
        - `['PF00102']`

    """

    def __init__(self, pfam_version='34.0', query=[['rcsb_polymer_entity_annotation.type','exact_match','Pfam']], **kwargs):
        self.pfam_version = pfam_version
        super().__init__(query=query, **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["protein"]["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        pfams = []
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'Pfam':
                # pfams.append(a['name'])
                pfams.append(a['annotation_id'])
        protein['protein']['Pfam'] = pfams
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Protein Family (Pfam)"
        desc['values'] = f"{len(set((p['Pfam'][0] for p in self.proteins)))} (root)"
        desc['type'] = 'Categorical, Hierarchical'
        return desc
