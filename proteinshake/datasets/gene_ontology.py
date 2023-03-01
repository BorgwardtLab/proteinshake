import json, os
from goatools.obo_parser import GODag

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import download_url

class GeneOntologyDataset(RCSBDataset):
    """ Proteins from RCSB for which the Gene Ontology (GO) term is known.
    Each protein in the dataset has a `GO` attribute which stores the path
    from the root to the leaves along the GO hierarchy.
    """

    additional_files = ['go-basic.obo']

    def __init__(self, query=[['rcsb_polymer_entity_annotation.type','exact_match','GO']], **kwargs):
        super().__init__(query=query, **kwargs)
        self.godag = GODag(f'{self.root}/go-basic.obo', prt=None)

    def download(self):
        super().download()
        if not os.path.exists(f'{self.root}/godag.obo'):
            download_url(f'http://current.geneontology.org/ontology/go-basic.obo', f'{self.root}', log=False)
        self.godag = GODag(f'{self.root}/go-basic.obo', prt=None)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["protein"]["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        go_terms = []
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'GO':
                go_terms.extend([go['id'] for go in a['annotation_lineage']])
        protein['protein']['molecular_function'] = [term for term in go_terms if self.godag[term].namespace == 'molecular_function']
        protein['protein']['cellular_component'] = [term for term in go_terms if self.godag[term].namespace == 'cellular_component']
        protein['protein']['biological_process'] = [term for term in go_terms if self.godag[term].namespace == 'biological_process']
        return protein

    def describe(self):
        desc = super().describe()
        desc['property'] = "Gene Ontology (GO)"
        desc['values'] = f"{len(set((p['GO'][0] for p in self.proteins)))} (root)"
        desc['type'] = 'Categorical, Hierarchical'
        return desc
