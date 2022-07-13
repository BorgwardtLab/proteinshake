import json

from torch_pdb.datasets import RCSBDataset

class ECDataset(RCSBDataset):
    """ Proteins from RCSB for which Enzyme Classification is known.
    Each item in this dataset has the attribute `EC` which is a string
    identifier.
    """

    def __init__(self, query=[['rcsb_polymer_entity.rcsb_ec_lineage.name','exists']], **kwargs):
        """

        Args:
            query: REST-API query.

        """
        super().__init__(query=query, **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        protein['EC'] = annot['rcsb_polymer_entity']['rcsb_ec_lineage'][-1]['id']
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = 'Enzyme Classification (`EC`)'
        desc['values'] = len(set((p['EC'] for p in self.proteins)))
        desc['type'] = 'Categorical'
        return desc
