import pandas as pd

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import download_url

class SCOPDataset(RCSBDataset):
    """ Proteins for which the SCOP classification is known.
    """

    def __init__(self, test=False, **kwargs):
        """

        Args:
            query: REST-API query.

        """
        query = [
                 ['rcsb_polymer_instance_annotation.type','exact_match', 'SCOP'],
                 ]
        self.test = test
        super().__init__(query=query, only_single_chain=True, **kwargs)


    def _parse_scop(self, path):
        names = ['FA-DOMID', 'FA-PDBID', 'FA-PDBREG', 'FA-UNIID', 'FA-UNIREG', 'SF-DOMID', 'SF-PDBID', 'SF-PDBREG', 'SF-UNIID', 'SF-UNIREG', 'SCOPCLA']
        return pd.read_csv(path, sep=' ', comment='#', names=names)

    def download_limit(self):
        if self.test:
            return 5
        else:
            return None

    def download(self):
        # get the proteins
        super().download()
        # get the annots
        download_url(f'http://scop.mrc-lmb.cam.ac.uk/files/scop-cla-latest.txt', f'{self.root}/raw/scop.txt')
        self.scop= self._parse_scop(f'{self.root}/raw/scop.txt')

    def add_protein_attributes(self, protein):
        """ We annotate the protein with the scop classifications at each level.

        SCOPCLA - SCOP domain classification. The abbreviations denote: TP=protein type, CL=protein class, CF=fold, SF=superfamily, FA=family
        """
        try:
            match = self.scop.loc[self.scop['FA-PDBID'] == protein['protein']['id']].iloc[0]
            print(match)
        except KeyError:
            return None
        scop_classes = dict([cla.split("=") for cla in match['SCOPCLA'].split(",")])
        for cla, val in scop_classes.keys():
            protein['protein'][cla] = val
        return protein
