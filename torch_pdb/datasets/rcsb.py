import requests, glob, json
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

from torch_pdb.datasets import TorchPDBDataset
from torch_pdb.utils import download_url

class RCSBDataset(TorchPDBDataset):
    """ Non-redundant structures taken from RCSB Protein Databank.
    """

    def __init__(self, query=[], similarity_cutoff=70, **kwargs):
        self.similarity_cutoff = similarity_cutoff
        self.query = query
        super().__init__(only_single_chain=True, **kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb.gz')

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        if self.n_jobs == 1:
            print('Warning: Downloading an RCSB dataset with use_precompute = False is very slow. Consider increasing n_jobs.')
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
                            "parameters": {"operator": "exact_match", "value": "Protein (only)", "attribute": "rcsb_entry_info.selected_polymer_entity_types"}
                        },
                        {
                            "type": "terminal",
                            "service": "text",
                            "parameters": {"attribute": "rcsb_entry_info.deposited_polymer_entity_instance_count", "operator": "equals", "value": 1}
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
            try:
                response_dict = json.loads(r.text)
                ids.extend([x['identifier'].split('_')[0] for x in response_dict['result_set']])
            except:
                print('An error occured when querying RCSB.')
                print(r.text)
                exit()
            if total is None:
                total = response_dict['group_by_count']
            i += batch_size
        ids = ids[:self.download_limit()] # for testing

        failed = Parallel(n_jobs=self.n_jobs)(delayed(self.download_from_rcsb)(id) for id in tqdm(ids, desc='Downloading PDBs'))
        failed = [f for f in failed if not f is True]
        if len(failed)>0:
            print(f'Failed to download {len(failed)} PDB files.')


    def download_from_rcsb(self, id):
        try:
            r = requests.get(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{id}/1')
            obj = json.loads(r.text)
            download_url(f'https://files.rcsb.org/download/{id}.pdb.gz', f'{self.root}/raw/files', log=False)
            with open(f'{self.root}/raw/files/{id}.annot.json', 'w') as file:
                json.dump(obj, file)
            return True
        except KeyboardInterrupt:
            exit()
        except:
            return id
