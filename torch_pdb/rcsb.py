import requests, glob, torch, json
import pandas as pd
from torch_geometric.data import download_url
from torch_pdb import TorchPDBDataset
from tqdm import tqdm
from joblib import Parallel, delayed

class RCSBDataset(TorchPDBDataset):

    def __init__(self, query=[], similarity_cutoff=70, **kwargs):
        self.similarity_cutoff = similarity_cutoff
        self.query = query
        super().__init__(**kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb.gz')

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        total = None
        i = 0
        batch_size = 5000
        ids = []
        while total is None or total > i:
            payload = {
                "query": {
                    "type": "group",
                    'logical_operator': 'and',
                    'nodes': [
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {"operator": "exact_match", "value": "Protein", "attribute": "entity_poly.rcsb_entity_polymer_type"}
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {"attribute": "rcsb_entry_info.entity_count", "operator": "equals", "value": 1}
                        },
                        *[
                            {
                                "type": "terminal",
                                "service": "text",
                                "parameters": {k:v for k,v in zip(['attribute','operator','value'], q)}
                            }
                        for q in self.query],
                    ],
                },
                "request_options": {
                    "group_by": {
                        "aggregation_method": "sequence_identity",
                        "similarity_cutoff": self.similarity_cutoff,
                    },
                    "group_by_return_type": "representatives",
                    "paginate": {"start": i, "rows": i+batch_size}
                },
                "return_type": "polymer_entity"
            }
            r = requests.get(f'https://search.rcsb.org/rcsbsearch/v2/query?json={json.dumps(payload)}')
            r = json.loads(r.text)
            ids.extend([x['identifier'].split('_')[0] for x in r['result_set']])
            if total is None:
                total = r['group_by_count']
            i += batch_size
            print(f'\rQuerying {min(i,total)} of {total}', end='')
        print()

        _ = Parallel(n_jobs=10)(delayed(self.download_from_rcsb)(id) for id in tqdm(ids, desc='Downloading PDBs'))

        self.download_complete()

    def download_from_rcsb(self, id):
        try:
            r = requests.get(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{id}/1')
            obj = json.loads(r.text)
            download_url(f'https://files.rcsb.org/download/{id}.pdb.gz', f'{self.root}/raw/files', log=False)
            with open(f'{self.root}/raw/files/{id}.annot.json', 'w') as file:
                json.dump(obj, file)
        except KeyboardInterrupt:
            exit()
        except:
            print(f'Downloading PDB ID {id} failed.')





class GODataset(RCSBDataset):

    def __init__(self, **kwargs):
        super().__init__(query=[['rcsb_polymer_entity_annotation.type','exact_match','GO']], **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        go_terms = []
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'GO':
                go_terms.extend([go['id'] for go in a['annotation_lineage']])
        protein['GO'] = go_terms
        return protein


class ECDataset(RCSBDataset):

    def __init__(self, **kwargs):
        super().__init__(query=[['rcsb_polymer_entity.rcsb_ec_lineage.name','exists']], **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        protein['EC'] = annot['rcsb_polymer_entity']['rcsb_ec_lineage'][-1]['id']
        return protein


class PfamDataset(RCSBDataset):

    def __init__(self, **kwargs):
        super().__init__(query=[['rcsb_polymer_entity_annotation.type','exact_match','Pfam']], **kwargs)

    def add_protein_attributes(self, protein):
        with open(f'{self.root}/raw/files/{protein["ID"]}.annot.json','r') as file:
            annot = json.load(file)
        pfams = []
        for a in annot['rcsb_polymer_entity_annotation']:
            if a['type'] == 'Pfam':
                pfams.append(a['name'])
        protein['Pfam'] = pfams
        return protein
