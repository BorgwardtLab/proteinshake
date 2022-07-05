from torch_pdb.datasets import RCSBDataset

class ECDataset(RCSBDataset):

    def __init__(self, query=[['rcsb_polymer_entity.rcsb_ec_lineage.name','exists']], **kwargs):
        super().__init__(query=query, **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        protein['EC'] = annot['rcsb_polymer_entity']['rcsb_ec_lineage'][-1]['id']
        return protein


