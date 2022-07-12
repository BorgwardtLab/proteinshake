import json

from torch_pdb.datasets import RCSBDataset

class GODataset(RCSBDataset):
    """ Proteins from RCSB for which the Gene Ontology (GO) term is known.
    Each protein in the dataset has a `GO` attribute which stores the path
    from the root to the leaves along the GO hierarchy.
    """

    def __init__(self, query=[['rcsb_polymer_entity_annotation.type','exact_match','GO']], **kwargs):
        super().__init__(query=query, **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        go_terms = []
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'GO':
                go_terms.extend([go['id'] for go in a['annotation_lineage']])
        protein['GO'] = go_terms
        return protein
    def describe(self):
        desc = super().describe()
        desc['property'] = "Gene Ontology (GO)"
        desc['values'] = f"{len(set(s[0] for s in self.data.GO))} (root)"
        desc['type'] = 'Categorical, Hierarchical'
        return desc

