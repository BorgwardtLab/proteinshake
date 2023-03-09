import json, os
from goatools.obo_parser import GODag

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import download_url, unzip_file
from functools import cached_property

class GeneOntologyDataset(RCSBDataset):
    """ Proteins from RCSB for which the Gene Ontology (GO) term is known.
    Each protein in the dataset has a `GO` attribute which stores the path
    from the root to the leaves along the GO hierarchy.
    """

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
        if not os.path.exists(f'{self.root}/godag.obo'):
            download_url(f'http://current.geneontology.org/ontology/go-basic.obo', f'{self.root}', log=False)
            os.rename(f'{self.root}/go-basic.obo', f'{self.root}/{self.name}.godag.obo')

    def add_protein_attributes(self, protein):
        godag = GODag(f'{self.root}/{self.name}.godag.obo', prt=None)
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
