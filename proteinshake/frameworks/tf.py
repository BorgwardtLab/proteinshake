import tensorflow as tf
from proteinshake.frameworks.dataset import FrameworkDataset

class TensorflowVoxelDataset(FrameworkDataset):

    def convert_to_framework(self, data_item):
        return tf.sparse.from_dense(tf.convert_to_tensor(data_item.data, dtype=tf.float32))

    def load_transform(self, data, protein_dict):
        return tf.sparse.to_dense(data), protein_dict


class TensorflowPointDataset(FrameworkDataset):

    def convert_to_framework(self, data_item):
        return tf.convert_to_tensor(data_item.data, dtype=tf.float32)
