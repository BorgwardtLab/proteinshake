# -*- coding: utf-8 -*-
import os, glob
import pandas as pd
import numpy as np
from joblib import Parallel, delayed
from biopandas.pdb import PandasPdb
from sklearn.neighbors import KDTree

from proteinshake.datasets import Dataset
from proteinshake.utils import extract_tar, download_url, progressbar, load, save, unzip_file

class ProteinProteinInterfaceDataset(Dataset):
    """Protein-protein complexes from PDBBind with annotated interfaces.
    Residues and atoms in each protein are marked with a boolean `is_interface` to indicate residues/atoms defined to belong to the interface of two protein chains.
    The default threshold for determining interface residues is 6 Angstroms (used by DIPS).
    See :meth:`proteinshake.utils.get_interfaces` for details.

    .. admonition:: Please cite

      Wang, Renxiao, et al. "The PDBbind database: Collection of binding affinities for protein− ligand complexes with known three-dimensional structures." Journal of medicinal chemistry 47.12 (2004): 2977-2980.

    .. admonition:: Source

      Raw data was obtained and modified with permission from `PDBbind-CN <http://www.pdbbind.org.cn/>`_, originally licensed under the `End User Agreement for Access to the PDBbind-CN Database and Web Site <http://www.pdbbind.org.cn/enroll.php>`_.


    Parameters
    ----------
    root: str
        Root directory where the dataset should be saved.
    name: str
        The name of the dataset.
    version: str
        PDBBind database version to use.
    cutoff: float
        Distance in angstroms within which a pair of residues is considered to
        belong to the interface.

   .. list-table:: Annotations
      :widths: 25 55 20
      :header-rows: 1

      * - Attribute
        - Key
        - Sample value
      * - Chain Identifier 
        - :code:`protein[{'residue' | 'atom'}]['chain_id']`
        - ``['X', 'X', ..'Y', 'Y']``
      * - Binding interface
        - :code:`protein[{'residue' | 'atom'}]['is_interface']`
        - ``[0, 0, .., 1, 0]``


    """
    additional_files = [
        'ProteinProteinInterfaceDataset.interfaces.json',
    ]

    def __init__(self, cutoff=6, version='2020', **kwargs):
        self.version = version
        self.cutoff = cutoff
        super().__init__(**kwargs)
        self.interfaces = load(f'{self.root}/{self.name}.interfaces.json')

    def get_raw_files(self):
        return sorted(glob.glob(f'{self.root}/raw/files/chains/*.pdb'))
    
    def get_id_from_filename(self, filename):
        return filename.rstrip('.pdb')#.rstrip('.ent')

    def get_contacts(self, path):
        """ Obtain interfacing residues within a single structure of polymers.
        Uses KDTree data structure for vector search.
        """
        df = PandasPdb().read_pdb(path).df['ATOM']
        
        # split chains
        pdbid = os.path.basename(path).rstrip('.pdb')
        for chain, chain_df  in df.groupby('chain_id'):
            new_df = PandasPdb()
            new_df._df = {'ATOM': chain_df}
            new_df.to_pdb(f'{self.root}/raw/files/chains/{pdbid}_{chain}.pdb')
            
        # compute contacts
        df = df[df['atom_name'] == 'CA']
        coords = df[['x_coord','y_coord','z_coord']].to_numpy()
        chain_ids = np.array(df['chain_id'])
        chain_index = np.hstack([np.arange(np.sum(chain_ids==chain)) for chain in list(dict.fromkeys(chain_ids))]) # get residue position relative to chain
        contacts = KDTree(coords, leaf_size=1).query_radius(coords, self.cutoff) # for each residue, get all neighbouring residues within the cutoff radius
        query_chains = np.repeat(chain_ids, [len(c) for c in contacts]) # construct long index of query chain ids
        query_chain_index = np.repeat(chain_index, [len(c) for c in contacts]) # construct long index of query chain index
        result_chains = np.hstack([chain_ids[c] for c in contacts]) # construct long index of result chain ids
        result_chain_index = np.hstack([chain_index[c] for c in contacts]) # construct long index of result chain index
        interfaces = pd.DataFrame({'query': query_chains, 'result': result_chains, 'index': tuple(zip(query_chain_index.tolist(),result_chain_index.tolist()))})
        interfaces = interfaces[query_chains != result_chains] # remove contacts in the same chain
        if len(interfaces) == 0: return pdbid, None
        interfaces = dict(interfaces.groupby('query').apply(lambda x: dict(x.groupby('result')['index'].apply(list)))) # format to dict of dicts: {chain_A: chain_B: [13, 27, ...]}
        return pdbid, interfaces

    def download(self):
        #download_url(f'https://pdbbind.oss-cn-hangzhou.aliyuncs.com/download/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw')
        #extract_tar(f'{self.root}/raw/PDBbind_v{self.version}_PP.tar.gz', f'{self.root}/raw/files', extract_members=True)
        os.makedirs(f'{self.root}/raw/files/chains', exist_ok=True)
        complexed_files = sorted(glob.glob(f'{self.root}/raw/files/PP/*.pdb'))
        contacts = Parallel(n_jobs=self.n_jobs)(delayed(self.get_contacts)(path) for path in progressbar(complexed_files, desc='Computing interfaces', verbosity=self.verbosity))
        interfaces = {pdbid: interfaces for pdbid, interfaces in contacts if not interfaces is None}
        save(interfaces, f'{self.root}/{self.name}.interfaces.json')
        
            
        
