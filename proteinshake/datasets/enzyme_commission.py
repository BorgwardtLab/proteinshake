import json

from proteinshake.datasets import RCSBDataset

class EnzymeCommissionDataset(RCSBDataset):
    """ Enzymes with annotated enzyme commission (EC) numbers.

    .. admonition:: Please cite

        Berman, H M et al. “The Protein Data Bank.” Nucleic acids research vol. 28,1 (2000): 235-42. doi:10.1093/nar/28.1.235

    .. admonition:: Source

        Raw data was obtained and modified from `RCSB Protein Data Bank <https://www.rcsb.org/>`_, originally licensed under `CC0 1.0 <https://creativecommons.org/publicdomain/zero/1.0/>`_.


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

    description = 'Enzymes'

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
