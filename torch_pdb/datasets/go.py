from torch_pdb.datasets import RCSBDataset

class GODataset(RCSBDataset):

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


