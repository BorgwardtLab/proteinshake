import os
import re
import tarfile
import glob

from tqdm import tqdm

from torch_pdb.datasets import TorchPDBDataset
from torch_pdb.utils import download_url, extract_tar, load, save

AF_DATASET_NAMES = {
    'arabidopsis_thaliana': 'UP000006548_3702_ARATH',
    'caenorhabditis_elegans': 'UP000001940_6239_CAEEL',
    'candida_albicans': 'UP000000559_237561_CANAL',
    'danio_rerio': 'UP000000437_7955_DANRE',
    'dictyostelium_discoideum': 'UP000002195_44689_DICDI',
    'drosophila_melanogaster': 'UP000000803_7227_DROME',
    'escherichia_coli': 'UP000000625_83333_ECOLI',
    'glycine_max': 'UP000008827_3847_SOYBN',
    'homo_sapiens': 'UP000005640_9606_HUMAN',
    'methanocaldococcus_jannaschii': 'UP000000805_243232_METJA',
    'mus_musculus': 'UP000000589_10090_MOUSE',
    'oryza_sativa': 'UP000059680_39947_ORYSJ',
    'rattus_norvegicus': 'UP000002494_10116_RAT',
    'saccharomyces_cerevisiae': 'UP000002311_559292_YEAST',
    'schizosaccharomyces_pombe': 'UP000002485_284812_SCHPO',
    'zea_mays': 'UP000007305_4577_MAIZE',
    'swissprot': 'swissprot_pdb',
}

class AlphaFoldDataset(TorchPDBDataset):
    """ 3D structures predicted from sequence by AlphaFold.
    """

    def __init__(self, organism, **kwargs):
        if organism == 'all':
            self.organism = [o for o in AF_DATASET_NAMES.keys() if o != 'swissprot']
        else:
            if type(organism) == str:
                self.organism = organism.lower().replace(' ','_')
            elif type(organism) == list:
                self.organism = [o.lower().replace(' ','_') for o in organism]
        self.base_url = 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/'
        super().__init__(**kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/*/*.pdb.gz')[:self.download_limit()]

    def get_id_from_filename(self, filename):
        return re.search('(?<=AF-)(.*)(?=-F1-model)', filename).group()

    def download_precomputed(self):
        # overload this to compile multiple organisms into one
        if os.path.exists(f'{self.root}/{self.__class__.__name__}.pt'):
            return
        def _download(organism):
            download_url(f'https://github.com/BorgwardtLab/torch-pdb/releases/download/{self.release}/{self.__class__.__name__}_{organism}.pt', f'{self.root}')
        if type(self.organism) == str:
            _download(self.organism)
            os.rename(f'{self.root}/{self.__class__.__name__}_{self.organism}.pt', f'{self.root}/{self.__class__.__name__}.pt')
        elif type(self.organism) == list:
            _ = [_download(organism) for organism in self.organism]
            proteins = [p for organism in self.organism for p in load(f'{self.root}/{self.__class__.__name__}_{organism}.pt')]
            save(proteins, f'{self.root}/{self.__class__.__name__}.h5')
            _ = [os.remove(f'{self.root}/{self.__class__.__name__}_{organism}.pt') for organism in self.organism]


    def download(self):
        def _download(organism):
            os.makedirs(f'{self.root}/raw/{organism}', exist_ok=True)
            download_url(self.base_url+AF_DATASET_NAMES[organism]+'_v2.tar', f'{self.root}/raw/{organism}')
            extract_tar(f'{self.root}/raw/{organism}/{AF_DATASET_NAMES[organism]}_v2.tar', f'{self.root}/raw/{organism}')
        if type(self.organism) == str:
            _download(self.organism)
        elif type(self.organism) == list:
            _ = [_download(organism) for organism in self.organism]
