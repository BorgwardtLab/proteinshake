import os
import re
import tarfile
import glob

from proteinshake.datasets import Dataset
from proteinshake.utils import download_url, extract_tar, load, save, unzip_file, progressbar

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

description = 'Predicted structures'

class AlphaFoldDataset(Dataset):
    """ SwissProt 3D structures predicted by AlphaFold.

    .. admonition:: Please cite

      Jumper, John, et al. "Highly accurate protein structure prediction with AlphaFold." Nature 596.7873 (2021): 583-589.

      Varadi, Mihaly, et al. "AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models." Nucleic acids research 50.D1 (2022): D439-D444.

    .. admonition:: Source

      Raw data was obtained and modified from `AlphaFoldDB <https://alphafold.ebi.ac.uk>`_, originally licensed under `CC-BY-4.0 <https://creativecommons.org/licenses/by/4.0/>`_.

    Parameters
    ----------
    organism: str
        The organism name or 'swissprot'.
    """

    exlude_args_from_signature = ['organism']

    def __init__(self, organism='swissprot', version='v4', only_single_chain=True, **kwargs):
        self.organism = organism.lower().replace(' ','_')
        self.base_url = 'https://ftp.ebi.ac.uk/pub/databases/alphafold/latest/'
        self.version = version
        super().__init__(only_single_chain=only_single_chain, **kwargs)

    @property
    def name(self):
        return f'{self.__class__.__name__}_{self.organism}'

    def get_raw_files(self):
        return glob.glob(f'{self.root}/raw/*/*.pdb')[:self.limit]

    def get_id_from_filename(self, filename):
        return re.search('(?<=AF-)(.*)(?=-F.+-model)', filename).group()

    def download(self):
        os.makedirs(f'{self.root}/raw/{self.organism}', exist_ok=True)
        download_url(self.base_url+AF_DATASET_NAMES[self.organism]+f'_{self.version}.tar', f'{self.root}/raw/{self.organism}', verbosity=self.verbosity)
        extract_tar(f'{self.root}/raw/{self.organism}/{AF_DATASET_NAMES[self.organism]}_{self.version}.tar', f'{self.root}/raw/{self.organism}', verbosity=self.verbosity)
        [unzip_file(f) for f in progressbar(glob.glob(f'{self.root}/raw/*/*.pdb.gz')[:self.limit], desc='Unzipping', verbosity=self.verbosity)]
