import pandas as pd
import torch
from biopandas.pdb import PandasPdb
from sklearn.neighbors import KDTree

import numpy as np

def get_interfaces(protein, ref_chain=None, cutoff=6):
    """Obtain interfacing residues within a single structure of polymers. Uses
    KDTree data structure for vector search. If structure not found in complex
    databse, it is automatically downloaded from RCSB.
    Parameters
    ----------
    protein: dict
        Parsed protein dictionary.
    Returns
    --------
        `list`: containing all Residue objects belonging to interface. As pairs
        of residues (res_1, res_2)
        `Structure`: BioPython Structure object
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
        resi_interface.append(r.item() in interface)
    return torch.tensor(resi_interface)

if __name__ == "__main__":
    df = PandasPdb().fetch_pdb('1pd7').df['ATOM']
    protein = {
        'sequence': ''.join(df['residue_name']),
        'residue_index': torch.tensor(df['residue_number'].tolist()).int(),
        'chain_id': df['chain_id'].tolist(),
        'coords': torch.tensor(df.apply(lambda row: (row['x_coord'], row['y_coord'], row['z_coord']), axis=1).to_list()).long(),
    }
    interface = get_interfaces(protein)
    print(interface)
    pass
