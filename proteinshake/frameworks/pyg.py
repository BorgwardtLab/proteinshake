import os
from itertools import repeat
import torch
from torch_geometric.utils import from_scipy_sparse_matrix
from torch_geometric.data import Data, InMemoryDataset
from proteinshake.utils import fx2str
from fastavro import reader as avro_reader
from tqdm import tqdm


def info2pyg(info):
    """ Try to convert info dictionary to tensor."""
    new_info = {}
    for k,v in info.items():
        if k in ['x','y','z']:
            continue
        try:
            new_info[k] = torch.Tensor(v).long()
        except:
            new_info[k] = v
            pass
    return new_info

class PygGraphDataset(InMemoryDataset):
    def __init__(self, graphDataset, transform=None, pre_transform=None, pre_filter=None):
        super().__init__(None, transform, pre_transform, pre_filter)
        path = f'{graphDataset.dataset.root}/processed/graph/{graphDataset.name}.pyg'
        transforms_repr = fx2str(transform)+fx2str(pre_transform)+fx2str(pre_filter)
        if not os.path.exists(path):
            data_list = []
            with open(f'{graphDataset.dataset.root}/{graphDataset.dataset.__class__.__name__}.{graphDataset.resolution}.avro', 'rb') as file:
                reader = avro_reader(file)
                for protein in tqdm(reader, total=int(reader.metadata['number_of_proteins']), desc='Converting to graph'):
                    nodes, adj = graphDataset.convert(protein)
                    nodes = torch.from_numpy(nodes)
                    edges = from_scipy_sparse_matrix(adj)
                    data = Data(x=nodes, edge_index=edges[0].long(), edge_attr=edges[1].unsqueeze(1).float(), **info2pyg(protein['protein']), **info2pyg(protein[graphDataset.resolution]))
                    data_list.append(data)
            if self.pre_filter is not None:
                data_list = [data for data in data_list if self.pre_filter(data)]
            if self.pre_transform is not None:
                data_list = [self.pre_transform(data) for data in data_list]
            data, slices = self.collate(data_list)
            torch.save((transforms_repr, data, slices), path)
        original_repr, self.data, self.slices = torch.load(path)
        assert original_repr == transforms_repr, f'The transforms are not the same as when the dataset was created. If you want to change them, delete the file at {path}'

    def get(self, idx):
        data = Data()
        for key in self.data.keys:
            item, slices = self.data[key], self.slices[key]
            if isinstance(item, torch.Tensor):
                s = list(repeat(slice(None), item.dim()))
                s[data.__cat_dim__(key, item)] = slice(slices[idx], slices[idx + 1])
                data[key] = item[s].clone()
            else:
                data[key] = [item[s] for s in range(slices[idx], slices[idx + 1])]
        return data
