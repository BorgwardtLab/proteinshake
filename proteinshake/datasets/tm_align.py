# -*- coding: utf-8 -*-
'''
TMalign needs to be in your $PATH. Follow the instructions at https://zhanggroup.org/TM-align/readme.c++.txt
'''
import glob
import requests
import os
import itertools
import re
import subprocess
import tempfile
import shutil
import numpy as np
from biopandas.pdb import PandasPdb
from collections import defaultdict
from joblib import Parallel, delayed
from tqdm import tqdm
from functools import cached_property

from proteinshake.datasets import RCSBDataset
from proteinshake.utils import (extract_tar,
                                download_url,
                                save,
                                load,
                                unzip_file,
                                global_distance_test,
                                local_distance_difference_test
                                )


class TMAlignDataset(RCSBDataset):
    """Proteins that were aligned with TMalign. The dataset provides the TM-score, RMSD, Global Distance Test (GDT), and Local Distance Difference Test (LDDT) as similarity/distance metrics between any two proteins.

    .. list-table:: Dataset stats
       :widths: 100
       :header-rows: 1

       * - # proteins
       * - 994


    .. code-block:: python

        from proteinshake.datasets import TMAlignDataset

        dataset = TMAlignDataset()
        proteins = dataset.proteins()
        protein_1, protein_2 = next(proteins)['protein']['ID'], next(proteins)['protein']['ID']

        dataset.tm_score(protein_1, protein_2)
        >>> 0.03
        dataset.rmsd(protein_1, protein_2)
        >>> 3.64
        dataset.gdt(protein_1, protein_2)
        >>> 0.61
        dataset.lddt(protein_1, protein_2)
        >>> 0.65

    """

    additional_files = [
        'TMAlignDataset.tmscore.npy',
        'TMAlignDataset.gdt.npy',
        'TMAlignDataset.rmsd.npy',
        'TMAlignDataset.lddt.npy'
    ]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if not self.use_precomputed: self.align_structures()
        self.protein_ids = [p['protein']['ID'] for p in self.proteins()]

        def download_file(filename):
            if not os.path.exists(f'{self.root}/{filename}'):
                download_url(f'{self.repository_url}/{filename}.gz', f'{self.root}', log=False)
                unzip_file(f'{self.root}/{filename}.gz')
            return load(f'{self.root}/{filename}')

        self._tm_score = download_file(f'{self.name}.tmscore.npy')
        self._rmsd = download_file(f'{self.name}.rmsd.npy')
        self._gdt = download_file(f'{self.name}.gdt.npy')
        self._lddt = download_file(f'{self.name}.lddt.npy')
        
    @property
    def limit(self):
        return 1000
    
    def align_structures(self):
        """ Calls TMAlignn on all pairs of structures and saves the output"""
        if os.path.exists(f'{self.root}/{self.name}.tmscore.npy'):
            return
        pdbids = [p['protein']['ID'] for p in self.proteins()]
        path_dict = {self.get_id_from_filename(os.path.basename(f)):f for f in self.get_raw_files()}
        paths = [path_dict[id] for id in pdbids]
        num_proteins = len(paths)
        combinations = np.array(list(itertools.combinations(range(num_proteins), 2)))
        TM, RMSD, GDT, LDDT = [np.ones((num_proteins,num_proteins), dtype=np.float16) * np.nan for _ in ['tm','rmsd','gdt','lddt']]
        np.fill_diagonal(TM, 1.0), np.fill_diagonal(RMSD, 0.0), np.fill_diagonal(GDT, 1.0), np.fill_diagonal(LDDT, 1.0)
        d = Parallel(n_jobs=self.n_jobs)(delayed(tmalign_wrapper)(paths[i], paths[j]) for i,j in tqdm(combinations, desc='Aligning'))
        x,y = tuple(combinations[:,0]), tuple(combinations[:,1])
        TM[x,y] = [x['TM1'] for x in d]
        TM[y,x] = [x['TM2'] for x in d]
        RMSD[x,y] = [x['RMSD'] for x in d]
        RMSD[y,x] = [x['RMSD'] for x in d]
        GDT[x,y] = [x['GDT'] for x in d]
        GDT[y,x] = [x['GDT'] for x in d]
        LDDT[x,y] = [x['LDDT'] for x in d]
        LDDT[y,x] = [x['LDDT'] for x in d]
        # save
        np.save(f'{self.root}/{self.name}.tmscore.npy', TM)
        np.save(f'{self.root}/{self.name}.rmsd.npy', RMSD)
        np.save(f'{self.root}/{self.name}.gdt.npy', GDT)
        np.save(f'{self.root}/{self.name}.lddt.npy', LDDT)

    def tm_score(self, protein_1, protein_2):
        return self._tm_score[self.protein_ids.index(protein_1)][self.protein_ids.index(protein_2)]
    
    def rmsd(self, protein_1, protein_2):
        return self._rmsd[self.protein_ids.index(protein_1)][self.protein_ids.index(protein_2)]

    def gdt(self, protein_1, protein_2):
        return self._gdt[self.protein_ids.index(protein_1)][self.protein_ids.index(protein_2)]

    def lddt(self, protein_1, protein_2):
        return self._lddt[self.protein_ids.index(protein_1)][self.protein_ids.index(protein_2)]

def tmalign_wrapper(pdb1, pdb2):
    """Compute TM score with TMalign between two PDB structures.
    Parameters
    ----------
    pdb1: str
        Path to PDB.
    pdb2 : str
        Path to PDB.
    return_superposition: bool
        If True, returns a protein dataframe with superposed structures.
    Returns
    -------
    dict
        Metric values TM1/TM2 (TM-Scores normalized to pdb1 or pdb2), RMSD, GDT
    """
    assert shutil.which('TMalign') is not None,\
           "No TMalign installation found. Go here to install : https://zhanggroup.org/TM-align/TMalign.cpp"
    with tempfile.TemporaryDirectory() as tmpdir:
        lines = subprocess.run(['TMalign','-outfmt','-1', pdb1, pdb2, '-o', f'{tmpdir}/superposition'], stdout=subprocess.PIPE).stdout.decode().split('\n')
        TM1 = lines[7].split()[1]
        TM2 = lines[8].split()[1]
        RMSD = lines[6].split()[4][:-1]
        seq1, ali, seq2 = lines[12], lines[13], lines[14]
        i, j, alignmentA, alignmentB = 0, 0, [], []
        for s1,a,s2 in zip(seq1,ali,seq2):
            if a != ' ': alignmentA.append(i)
            if a != ' ': alignmentB.append(j)
            if s1 != '-': i += 1
            if s2 != '-': j += 1
        os.rename(f'{tmpdir}/superposition_all', f'{tmpdir}/superposition_all.pdb')
        superposition = PandasPdb().read_pdb(f'{tmpdir}/superposition_all.pdb')
        df = superposition.df['ATOM']
        A = df[df['chain_id'] == 'A']
        B = df[df['chain_id'] == 'B']
        coordsA = np.array(list(zip(A['x_coord'], A['y_coord'], A['z_coord'])))[alignmentA]
        coordsB = np.array(list(zip(B['x_coord'], B['y_coord'], B['z_coord'])))[alignmentB]
        GDT = global_distance_test(coordsA, coordsB)
        LDDT = local_distance_difference_test(coordsA, coordsB)
    return {
        'TM1': float(TM1),
        'TM2': float(TM2),
        'RMSD': float(RMSD),
        'GDT': GDT,
        'LDDT': LDDT
    }
