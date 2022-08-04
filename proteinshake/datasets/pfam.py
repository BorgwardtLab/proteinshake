import json

from proteinshake.datasets import RCSBDataset

class PfamDataset(RCSBDataset):
    """ Proteins from RCSB for which the protein family (Pfam) term is known.
    Each protein in the dataset has a `Pfam` attribute which stores the list of protein families.
    """

    def __init__(self, query=[['rcsb_polymer_entity_annotation.type','exact_match','Pfam']], **kwargs):
        super().__init__(query=query, **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        pfams = []
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'Pfam':
                pfams.append(a['name'])
        protein['Pfam'] = pfams
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Protein Family (Pfam)"
        desc['values'] = f"{len(set((p['Pfam'][0] for p in self.proteins)))} (root)"
        desc['type'] = 'Categorical, Hierarchical'
        return desc
