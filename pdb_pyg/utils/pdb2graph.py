# -*- coding: utf-8 -*-
import itertools
import pandas as pd
import networkx as nx
import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import pdist, squareform
from sklearn.neighbors import kneighbors_graph, radius_neighbors_graph

from pyg_pdb.utils.resi_atoms import RESI_THREE_TO_1, AA_TO_INDEX
from pyg_pdb.utils.aaindex import get_aaindex_features


BOND_TYPES = {'backbone': 0, 'structure': 1}

def pdb2df(filename):
    atomic_df = PandasPdb().read_pdb(filename).df['ATOM']
    # sort the df
    atomic_df.sort_values('atom_number', inplace=True)
    return atomic_df


class ProteinGraph(object):
    """Parse a pdb/mmcif file to a graph
    """
    def __init__(self,
            datapath,
            uniprot_id=None,
            node_type="residue",
            neighbor_type="knn",
            knn=5,
            radius=5,
        ):
        self.datapath = datapath
        self.uniprot_id = uniprot_id
        self.node_type = node_type
        self.neighbor_type = neighbor_type
        self.knn = knn
        self.radius = radius

        self.atomic_df = pdb2df(datapath)

        self._initialize_graph()
        self._add_nodes_to_graph()
        self._add_edges_to_graph()
        self._add_aaindex_node_features()

    def get_sequences(self):
        """Output raw sequences from the atomic dataframe.
        """
        # TODO: can we accelerate it?
        # TODO: convert three to one!!!
        atomic_df = self.atomic_df
        sequences = atomic_df.loc[atomic_df['atom_name'] == 'CA']
        # sequences = sequences.sort_values('residue_number')
        sequences = sequences.groupby(['chain_id'])['residue_name']
        sequences = sequences.apply(
            lambda x: x.apply(lambda x: RESI_THREE_TO_1[x]).str.cat()).tolist()
        return sequences

    def _initialize_graph(self):
        """Initialize the protein graph with meta data.
        """
        atomic_df = self.atomic_df
        chain_ids = list(atomic_df["chain_id"].unique())
        graph = nx.Graph(
            name=self.uniprot_id,
            chain_ids=chain_ids,
            raw_pdb=atomic_df,
        )
        self.graph = graph

    def _add_nodes_to_graph(self):
        """add nodes after initializing the graph
        """
        self.atomic_df["node_id"] = (
            self.atomic_df["chain_id"]
            + ":"
            + self.atomic_df["residue_number"].apply(str)
            + ":"
            + self.atomic_df["residue_name"]
        )
        if self.node_type == "residue":
            atomic_df = self.atomic_df.loc[self.atomic_df['atom_name'] == 'CA']
            self.graph.coords = atomic_df[['x_coord', 'y_coord', 'z_coord']].to_numpy()
            for residue_position, (node_id, chain_id, residue_name, residue_number, coord) in enumerate(zip(
                atomic_df['node_id'],
                atomic_df['chain_id'],
                atomic_df['residue_name'],
                atomic_df['residue_number'],
                atomic_df[['x_coord', 'y_coord', 'z_coord']].to_numpy(),
            )):
                self.graph.add_node(
                    node_id,
                    chain_id=chain_id,
                    residue_idx=AA_TO_INDEX[RESI_THREE_TO_1[residue_name]],
                    residue=RESI_THREE_TO_1[residue_name],
                    residue_name=residue_name,
                    residue_number=residue_number,
                    residue_position=residue_position,
                    # distance=dist,
                    # nn=nn,
                    coord=coord,
                )
        else:
            raise ValueError("Not implemented!")

    def _add_edges_to_graph(self):
        # add backbone edges
        df = pd.DataFrame(
            self.graph.nodes(data='chain_id'), columns=['node_id', 'chain_id'])
        for chain_id, node_id in df.groupby(['chain_id'])['node_id']:
            node_i, node_j = itertools.tee(node_id)
            next(node_j, None)
            # bond_type: backbone=0, structure=1
            self.graph.add_edges_from(zip(node_i, node_j), bond_type=BOND_TYPES['backbone'])

        # add distance-based edges
        # coords = self.graph.coords
        # node_ids = df['node_id']
        # if self.neighbor_type == "knn":
        #     dist = kneighbors_graph(coords, n_neighbors=self.knn, mode='connectivity')
        # elif self.neighbor_type == "radius":
        #     dist = radius_neighbors_graph(coords, radius=self.radius, mode='connectivity')

        # dist = dist.tocoo()
        # # print(dist)
        # # dist = squareform(pdist(coords))
        # for i, j in zip(dist.row, dist.col):
        #     node_i, node_j = node_ids[i], node_ids[j]
        #     # if self.graph.has_edge(node_i, node_j):
        #     #     self.graph.edges[node_i, node_j]['kind'].add("knn")
        #     if not self.graph.has_edge(node_i, node_j):
        #         self.graph.add_edge(node_i, node_j, bond_type=BOND_TYPES['structure'])


        # other types of edges

    def _add_aaindex_node_features(self):
        # todo: check for natural amino acids
        # todo: PCA on aaindex features
        # todo: impute missing values
        attrs = {node: {'aa_idx':get_aaindex_features(
            data['residue'])} for node, data in self.graph.nodes(data=True)}
        nx.set_node_attributes(self.graph, attrs)

    def get_graph(self):
        return self.graph


if __name__ == "__main__":
    pro = ProteinGraph("../datasets/AF-Q8W3K0-F1-model_v1.pdb",
            neighbor_type="radius")
    print(pro.get_sequences())
    # print(pro.get_graph())
    # print(pro.get_graph().edges(data='kind'))
    # print(pro.get_sequences())
