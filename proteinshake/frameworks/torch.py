import torch
from torch.utils.data import Dataset as TorchDataset
from proteinshake.frameworks.dataset import FrameworkDataset


class TorchVoxelDataset(FrameworkDataset, TorchDataset):

    def convert_to_framework(self, data_item):
        return torch.tensor(data_item.data).float().to_sparse()

    def load_transform(self, data, protein_dict):
        return data.to_dense(), protein_dict


class TorchPointDataset(FrameworkDataset, TorchDataset):

    def convert_to_framework(self, data_item):
        return torch.tensor(data_item.data).float()
