import requests, re, io, time
from requests.adapters import HTTPAdapter, Retry
import pandas as pd
import numpy as np
from proteinshake.utils import progressbar

def uniprot_query(query, columns='', verbosity=2):
    columns = 'accession,'+columns
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    def get_next_link(headers):
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)
    def get_batch(batch_url):
        while batch_url:
            response = session.get(batch_url)
            response.raise_for_status()
            total = response.headers["x-total-results"]
            yield response, total
            batch_url = get_next_link(response.headers)
    url = f'https://rest.uniprot.org/uniprotkb/search?fields={columns}&format=tsv&query={query}&size=500'
    df = pd.DataFrame()
    with progressbar(total=100, verbosity=verbosity) as pbar:
        for batch, total in get_batch(url):
            pbar.total = int(total)
            batch = pd.read_csv(io.StringIO(batch.text), sep='\t')
            df = pd.concat([df,batch])
            pbar.update(len(batch))
    df = df.set_index('Entry', drop=True).replace(np.nan, None)
    return df.to_dict('index')


def uniprot_map(ids, source, target, polling_interval=3):
    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))
    response = requests.post(f'https://rest.uniprot.org/idmapping/run',
        data={'from': source, 'to': target, 'ids': ','.join(ids)},
    )
    response.raise_for_status()
    job_id = response.json()['jobId']
    while True:
        response = session.get(f'https://rest.uniprot.org/idmapping/status/{job_id}')
        response.raise_for_status()
        response = response.json()
        if 'jobStatus' in response:
            if response['jobStatus'] == 'RUNNING': time.sleep(polling_interval)
            else: raise Exception(response['jobStatus'])
        elif bool(response['results'] or response['failedIds']):
            response = session.get(f'https://rest.uniprot.org/idmapping/stream/{job_id}')
            response.raise_for_status()
            results = response.json()['results']
            results = {r['from']:r['to'] for r in results[::-1]}
            mapped_ids = [results[id] if id in results else None for id in ids]
            return mapped_ids

