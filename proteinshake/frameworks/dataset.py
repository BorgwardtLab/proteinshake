import os
from proteinshake.utils import save, load, fx2str, progressbar, error

class FrameworkDataset():
    """ Dataset base class for different frameworks.

    Parameters
    ----------
    data_list: generator
        A generator of objects from a representation.
    size: int
        The size of the dataset.
    path: str
        Path to save the processed dataset.
    transform: function
        A transform function to be applied in the __getitem__ method. Signature: transform(data, protein_dict) -> (data, protein_dict)
    pre_transform: function
        A transform function to be applied before writing the data. Signature: transform((data, protein_dict)) -> (data, protein_dict)
    pre_filter: function
        A filter function to be applied before writing the data. Signature: transform(data, protein_dict) -> bool
    """

    def __init__(self, data_list, size, path, transform=None, pre_transform=None, pre_filter=None, verbosity=2):
        os.makedirs(path, exist_ok=True)
        self.verbosity = verbosity
        self.path = path
        self.transform = transform
        self.pre_transform = pre_transform
        self.pre_filter = pre_filter
        transforms_repr = fx2str(pre_transform) + fx2str(pre_filter)
        if not os.path.exists(f'{path}/{size-1}.pkl'):
            i = 0
            for data_item in progressbar(data_list, desc='Converting', total=size, verbosity=self.verbosity):
                data = self.convert_to_framework(data_item)
                protein_dict = data_item.protein_dict
                if not self.pre_filter is None and not self.pre_filter(data, protein_dict):
                    continue
                if not self.pre_transform is None:
                    data, protein_dict = self.pre_transform(data, protein_dict)
                save((data, protein_dict), f'{path}/{i}.pkl')
                i += 1
            save(i,f'{path}/size.pkl')
            save(transforms_repr,f'{path}/transforms.pkl')
        self.size = load(f'{path}/size.pkl')
        original_repr = load(f'{path}/transforms.pkl')
        if not original_repr == transforms_repr: error(f'The pre_transform and/or pre_filter are not the same as when the dataset was created. If you want to change them, delete the folder at {path}', verbosity=self.verbosity)

    def convert_to_framework(self, data_item):
        """ Converts data_item to a data object of the framework.
        """
        return data_item.data

    def load_transform(self, data, protein_dict):
        """ Applies a transform after loading, for example if the data has been stored in sparse format and needs to be converted to dense.
        """
        return data, protein_dict

    def __len__(self):
        return self.size

    def __getitem__(self, idx):
        try:
            idx = int(idx)
        except:
            return [self.__getitem__(i) for i in idx]
        if idx > self.size - 1:
            raise StopIteration
        data, protein_dict = self.load_transform(*load(f'{self.path}/{idx}.pkl'))
        if not self.transform is None:
            return self.transform((data, protein_dict))
        else:
            return data, protein_dict

    def len(self):
        pass

    def get(self):
        pass
