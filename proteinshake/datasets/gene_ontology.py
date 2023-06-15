import json, os
from goatools.obo_parser import GODag

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import download_url, unzip_file
from functools import cached_property

class GeneOntologyDataset(RCSBDataset):
    """ Proteins with annotated Gene Ontology (GO) terms.
    Each protein in the dataset has a `GO` attribute which stores the path from the root to the leaves along the GO hierarchy.
    The GeneOntologyDataset also has a `godag` attribute, which stores the GO hierarchy (see [goatools.obo_parser.GODag](https://github.com/tanghaibao/goatools)).

    .. admonition:: Please cite

      Botstein, David, et al. "Gene Ontology: tool for the unification of biology." Nat genet 25.1 (2000): 25-9.

      Central, G. O., et al. "The Gene Ontology knowledgebase in 2023." Genetics 224.1 (2023).

      Berman, H M et al. “The Protein Data Bank.” Nucleic acids research vol. 28,1 (2000): 235-42. doi:10.1093/nar/28.1.235

    .. admonition:: Source

      Raw data was obtained and modified from `RCSB Protein Data Bank <https://www.rcsb.org/>`_, originally licensed under `CC0 1.0 <https://creativecommons.org/publicdomain/zero/1.0/>`_.



    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 32633


   .. list-table:: Annotations
      :widths: 25 25 50
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Molecular Function
        - :code:`protein['protein']['molecular_function']`
        - :code:`['GO:0003674', 'GO:0005198']`
      * - Localization
        - :code:`protein['protein']['cellular_component']`
        - :code:`['GO:0005575', 'GO:0018995',..]`
      * - Biological process
        - :code:`protein['protein']['biological_process']`
        - ...


    """

    description = 'Gene Ontology'

    additional_files = ['GeneOntologyDataset.godag.obo']

    def __init__(self, query=[['rcsb_polymer_entity_annotation.type','exact_match','GO']], **kwargs):
        super().__init__(query=query, **kwargs)

    @cached_property
    def godag(self):
        if not os.path.exists(f'{self.root}/{self.name}.godag.obo'):
            download_url(f'{self.repository_url}/{self.name}.godag.obo.gz', f'{self.root}')
            unzip_file(f'{self.root}/{self.name}.godag.obo.gz')
        return GODag(f'{self.root}/{self.name}.godag.obo', prt=None)

    def download(self):
        super().download()
        if not os.path.exists(f'{self.root}/{self.name}.godag.obo'):
            download_url(f'http://current.geneontology.org/ontology/go-basic.obo', f'{self.root}', verbosity=0)
            os.rename(f'{self.root}/go-basic.obo', f'{self.root}/{self.name}.godag.obo')

    def add_protein_attributes(self, protein):
        godag = GODag(f'{self.root}/{self.name}.godag.obo', prt=None) # cannot use self.godag because the GODAG is not pickleable (for the release)
        with open(f'{self.root}/raw/files/{protein["protein"]["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        go_terms = []
        if not 'rcsb_polymer_entity_annotation' in annot: return None
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'GO':
                go_terms.extend([go['id'] for go in a['annotation_lineage']])
        protein['protein']['molecular_function'] = [term for term in go_terms if godag[term].namespace == 'molecular_function']
        protein['protein']['cellular_component'] = [term for term in go_terms if godag[term].namespace == 'cellular_component']
        protein['protein']['biological_process'] = [term for term in go_terms if godag[term].namespace == 'biological_process']
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Gene Ontology (GO)"
        desc['values'] = f"{len(set((p['GO'][0] for p in self.proteins)))} (root)"
        desc['type'] = 'Categorical, Hierarchical'
        return desc
