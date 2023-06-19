import glob

from proteinshake.datasets import Dataset
from proteinshake.utils import download_url, extract_tar, unzip_file, progressbar

class AlphaFoldDataset(Dataset):
    """ SwissProt 3D structures predicted by AlphaFold.

    .. admonition:: Please cite

      Jumper, John, et al. "Highly accurate protein structure prediction with AlphaFold." Nature 596.7873 (2021): 583-589.

      Varadi, Mihaly, et al. "AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models." Nucleic acids research 50.D1 (2022): D439-D444.

    .. admonition:: Source

      Raw data was obtained and modified from `AlphaFoldDB <https://alphafold.ebi.ac.uk>`_, originally licensed under `CC-BY-4.0 <https://creativecommons.org/licenses/by/4.0/>`_.

    Parameters
    ----------
    version: int, default 4
        The AlphaFoldDB version.
    """

    def __init__(self, version=4, only_single_chain=True, **kwargs):
        self.base_url = f'https://ftp.ebi.ac.uk/pub/databases/alphafold/v{version}'
        #self.file_name = f'swissprot_pdb_v{version}'
        self.file_name = f'UP000000805_243232_METJA_v{version}'
        super().__init__(only_single_chain=only_single_chain, **kwargs)

    def get_id_from_filename(self, filename):
        return filename.split('-')[1]

    def download(self):
        download_url(f'{self.base_url}/{self.file_name}.tar', f'{self.root}/raw', verbosity=self.verbosity)
        extract_tar(f'{self.root}/raw/{self.file_name}.tar', f'{self.root}/raw/files', verbosity=self.verbosity)
        for path in progressbar(glob.glob(f'{self.root}/raw/*/*.pdb.gz')[:self.limit], desc='Unzipping', verbosity=self.verbosity): unzip_file(path)
