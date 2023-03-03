import requests, glob, json, os, random
import pandas as pd
from tqdm import tqdm
from joblib import Parallel, delayed

from proteinshake.datasets import Dataset
from proteinshake.utils import download_url, unzip_file

class RCSBDataset(Dataset):
    """ Non-redundant structures taken from RCSB Protein Databank.

    This class also serves as a base class for all RCSB derived datasets. It can be subclassed by defining a default ``query`` argument. The query is a list of triplets ``(attribute, operator, value)`` according to `this <https://search.rcsb.org/#attribute-queries>`_ and `this <https://data.rcsb.org/data-attributes.html>`_ , which is passed to the REST API call to RCSB. See e.g. the GODataset subclass for an example. To find the right attributes, the queries can be constructed by doing an advanced search `at RCSB <https://www.rcsb.org/search/advanced>`_ and exporting to JSON. Also compare the API call in the :meth:`download()` method.

    It uses RCSB's integrated sequence similarity filtering to remove redundant proteins.

    Also, only single chain proteins are used. Change the REST payload in `download` to override this behaviour.

    Parameters
    ----------
    query: list
        A list of triplets `(attribute, operator, value)` to be added to the REST API call to RCSB.
    """

    def __init__(self, query=[], from_list=None, only_single_chain=True, max_requests=20, **kwargs):
        self.query = query
        self.from_list = from_list
        self.max_requests = max_requests
        super().__init__(only_single_chain=only_single_chain, **kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/files/*.pdb')

    def get_id_from_filename(self, filename):
        return filename[:4]

    def download(self):
        """ Fetches PDBs from RCSB with an API call. The default query selects protein-only
        structures with a single chain.
        """

        if self.n_jobs == 1:
            print('Warning: Downloading an RCSB dataset with use_precompute = False is very slow. Consider increasing n_jobs.')
        total = None
        i = 0
        batch_size = 5000
        ids = []
        if not self.from_list is None:
            ids = self.from_list
        else:
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
                            "similarity_cutoff": 100,
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

        ids = sorted(list(set(ids))) # filter identical ids
        random.seed(42)
        random.shuffle(ids) # for reproducible subsampling when using self.limit
        ids = ids[:self.limit] # for testing

        n_jobs = min(self.n_jobs, self.max_requests) # RCSB has a request limit
        if n_jobs < 1:
            n_jobs = self.max_requests

        failed = Parallel(n_jobs=n_jobs)(delayed(self.download_from_rcsb)(id) for id in tqdm(ids, desc='Downloading PDBs'))
        failed = [f for f in failed if not f is True]
        if len(failed)>0:
            print(f'Failed to download {len(failed)} PDB files.')


    def download_from_rcsb(self, id):
        try:
            r = requests.get(f'https://data.rcsb.org/rest/v1/core/polymer_entity/{id}/1')
            obj = json.loads(r.text)
            download_url(f'https://files.rcsb.org/download/{id}.pdb.gz', f'{self.root}/raw/files', log=False)
            unzip_file(f'{self.root}/raw/files/{id}.pdb.gz')
            with open(f'{self.root}/raw/files/{id}.annot.json', 'w') as file:
                json.dump(obj, file)
            return True
        except KeyboardInterrupt:
            exit()
        except Exception as e:
            if os.path.exists(f'{self.root}/raw/files/{id}.pdb.gz'):
                os.remove(f'{self.root}/raw/files/{id}.pdb.gz')
            return id
