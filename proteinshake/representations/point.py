import os
from tqdm import tqdm
from proteinshake.utils import checkpoint, one_hot, compose_embeddings

class PointDataset():
    """ Point cloud representation of a protein structure dataset.

    Converts a protein object to a point cloud.

    Parameters
    ----------
    embedding: Union[function, list]
        A function or list of functions for embedding the protein sequence to node attributes.

    """

    def __init__(self, root, proteins, embedding=one_hot):
        self.root = root
        if type(embedding) == list:
            self.embedding = compose_embeddings(embedding)
            emb_names = '_'.join([e.__name__ for e in embedding])
            self.name = f'emb_{emb_names}'
        else:
            self.embedding = embedding
            self.name = f'emb_{embedding.__name__}'
        self.proteins = self.convert(proteins)

    def protein2point(self, protein):
        protein['x'] = self.embedding(protein['sequence'])
        return protein

    @checkpoint('{root}/processed/point/{name}.pkl')
    def convert(self, proteins):
        return [self.protein2point(p) for p in tqdm(proteins, desc='Converting proteins to point clouds')]

    def torch(self):
        import torch
        from torch.utils.data import Dataset as TorchDataset
        class Dataset(TorchDataset):
            def __init__(self, path, proteins):
                if os.path.exists(path):
                    self.proteins = torch.load(path)
                else:
                    self.proteins = []
                    for p in proteins:
                        for key, value in p.items():
                            try:
                                p[key] = torch.tensor(value)
                            except:
                                pass
                        self.proteins.append(p)
                    torch.save(self.proteins, path)
            def __len__(self):
                return len(self.proteins)
            def __getitem__(self, i):
                return self.proteins[i]
        ds = Dataset(f'{self.root}/processed/point/{self.name}.torch.pkl', self.proteins)
        return ds

    def tf(self):
        import tensorflow as tf
        def generator():
            for p in self.proteins:
                yield p['coords']
        ds = tf.data.Dataset.from_generator(generator, tf.float32, output_shapes=[None,3])
        return ds
