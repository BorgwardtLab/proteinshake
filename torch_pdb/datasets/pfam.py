import json

from torch_pdb.datasets import RCSBDataset

class PfamDataset(RCSBDataset):

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

    def describe(self):
        desc = super().describe()
        desc['property'] = "Protein Family (Pfam)"
        desc['values'] = f"{len(set(s[0] for s in self.data.Pfam))} (root)"
        desc['type'] = 'Categorical, Hierarchical'
