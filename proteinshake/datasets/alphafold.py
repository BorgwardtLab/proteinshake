import os
import re
import tarfile
import glob

from tqdm import tqdm

from proteinshake.datasets import TorchPDBDataset
from proteinshake.utils import download_url, extract_tar, load, save, unzip_file

# A map of organism names to their download file names. See https://alphafold.ebi.ac.uk/download
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

    Requires the `organism` name to be specified. Can be a single organism, or a list of organism names, in which case the data will be concatenated. See https://alphafold.ebi.ac.uk/download for a full list of available organsims.
    Pass the full latin organism name separated by a space. `organism` can also be 'all', in which case the data of all organisms will be concatenated, or 'swissprot', in which case the full SwissProt structure predictions will be downloaded (ca. 500.000).

    Parameters
    ----------
    organism: Union[str, list]
        The organism name or a list of names or 'all' or 'swissprot'.
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
        super().__init__(only_single_chain=True, **kwargs)

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/*/*.pdb.gz')[:self.download_limit()]

    def get_id_from_filename(self, filename):
        return re.search('(?<=AF-)(.*)(?=-F.+-model)', filename).group()

    def download_precomputed(self, resolution='residue'):
        # overloads the parent method to compile multiple organisms into one
        if os.path.exists(f'{self.root}/{self.__class__.__name__}.residue.avro'):
            return
        def _download(organism):
            download_url(f'https://github.com/BorgwardtLab/torch-pdb/releases/download/{self.release}/{self.__class__.__name__}_{organism}.json.gz', f'{self.root}')
            print('Unzipping...')
            unzip_file(f'{self.root}/{self.__class__.__name__}_{organism}.json.gz')
        if type(self.organism) == str:
            _download(self.organism)
            os.rename(f'{self.root}/{self.__class__.__name__}_{self.organism}.json', f'{self.root}/{self.__class__.__name__}.json')
        elif type(self.organism) == list:
            _ = [_download(organism) for organism in self.organism]
            proteins = [p for organism in self.organism for p in load(f'{self.root}/{self.__class__.__name__}_{organism}.json')]
            save(proteins, f'{self.root}/{self.__class__.__name__}.json')
            _ = [os.remove(f'{self.root}/{self.__class__.__name__}_{organism}.json') for organism in self.organism]


    def download(self):
        def _download(organism):
            os.makedirs(f'{self.root}/raw/{organism}', exist_ok=True)
            download_url(self.base_url+AF_DATASET_NAMES[organism]+'_v3.tar', f'{self.root}/raw/{organism}')
            extract_tar(f'{self.root}/raw/{organism}/{AF_DATASET_NAMES[organism]}_v3.tar', f'{self.root}/raw/{organism}')
        if type(self.organism) == str:
            _download(self.organism)
        elif type(self.organism) == list:
            _ = [_download(organism) for organism in self.organism]
