import pandas as pd
from joblib import Parallel, delayed

from tqdm import tqdm

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import download_url

class SCOPDataset(RCSBDataset):
    """ Proteins for which the SCOP classification is known.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _parse_scop(self, path):
        names = ['FA-DOMID', 'FA-PDBID', 'FA-PDBREG', 'FA-UNIID', 'FA-UNIREG', 'SF-DOMID', 'SF-PDBID', 'SF-PDBREG', 'SF-UNIID', 'SF-UNIREG', 'SCOPCLA']
        return pd.read_csv(path, sep=' ', comment='#', names=names)

    def download(self):
        # get the annots
        download_url(f'http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt', f'{self.root}/raw/scop.txt')
        self.scop = self._parse_scop(f'{self.root}/raw/scop.txt')
        ids = list(self.scop['FA-PDBID'].unique())

        # get the proteins
        if self.n_jobs == 1:
            print('Warning: Downloading an RCSB dataset with use_precompute = False is very slow. Consider increasing n_jobs.')
        ids = ids[:self.limit] # for testing

        failed = Parallel(n_jobs=self.n_jobs)(delayed(self.download_from_rcsb)(id) for id in tqdm(ids, desc='Downloading PDBs'))
        failed = [f for f in failed if not f is True]
        if len(failed)>0:
            print(f'Failed to download {len(failed)} PDB files.')

    def add_protein_attributes(self, protein):
        """ We annotate the protein with the scop classifications at each level.

        SCOPCLA - SCOP domain classification. The abbreviations denote: TP=protein type, CL=protein class, CF=fold, SF=superfamily, FA=family
        """
        try:
            match = self.scop.loc[self.scop['FA-PDBID'] == protein['protein']['ID'].upper()].iloc[0]
        except KeyError:
            return None
        scop_classes = dict([cla.split("=") for cla in match['SCOPCLA'].split(",")])
        for cla, val in scop_classes.items():
            protein['protein']['SCOP-' + cla] = val
        return protein
