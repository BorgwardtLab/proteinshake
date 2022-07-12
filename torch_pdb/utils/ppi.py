import pandas as pd
from biopandas.pdb import PandasPdb
from sklearn.neighbors import KDTree

import numpy as np

def get_interfaces(protein, cutoff=6):
    """Obtain interfacing residues within a single structure of polymers. Uses
    KDTree data structure for vector search.

    Parameters
    ----------
    protein: dict
        Parsed protein dictionary.

    Returns
    --------
        `list`: indicator list for each residue with 0 if not in interface and 1 else.
    """

    #3-D KD tree
    df = pd.DataFrame({'x': [p[0] for p in protein['coords']],
                       'y': [p[1] for p in protein['coords']],
                       'z': [p[2] for p in protein['coords']],
                       'chain': protein['chain_id'],
                       'residue_index':protein['residue_index']
                       })
    resi_df = df.groupby(['residue_index', 'chain']).mean().reset_index()
    resi_coords = np.array([resi_df['x'].tolist(), resi_df['y'].tolist(), resi_df['z'].tolist()]).T
    kdt = KDTree(resi_coords, leaf_size=1)

    query = kdt.query_radius(resi_coords, cutoff)
    interface = set()
    for i,result in enumerate(query):
        res_index = resi_df.iloc[i].name
        this_chain = resi_df.iloc[i].chain
        for r in result:
            that_resi = resi_df.iloc[r].name
            that_chain = resi_df.iloc[r].chain
            if this_chain != that_chain:
                interface.add(res_index)
                interface.add(that_resi)

    resi_interface = []
    for r in protein['residue_index']:
        resi_interface.append(r in interface)
    return resi_interface
