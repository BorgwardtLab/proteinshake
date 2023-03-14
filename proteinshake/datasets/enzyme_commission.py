import json

from proteinshake.datasets import RCSBDataset

class EnzymeCommissionDataset(RCSBDataset):
    """
    Datsaet of enzymes which are assigned label based
    on the reaction they catalyze.


    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 15603 


   .. list-table:: Annotations
      :widths: 25 35 45
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Enzyme Commission
        - :code:`protein['protein']['EC']`
        - :code:`'2.7.7.4'`

    """

    def __init__(self, query=[['rcsb_polymer_entity.rcsb_ec_lineage.name','exists']], **kwargs):
        """

        Args:
            query: REST-API query.

        """
        super().__init__(query=query, **kwargs)

    def add_protein_attributes(self, protein):
        """ Fetch the enzyme class for each protein.
        """
        with open(f'{self.root}/raw/files/{protein["protein"]["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        protein['protein']['EC'] = annot['rcsb_polymer_entity']['rcsb_ec_lineage'][-1]['id']
        return protein
