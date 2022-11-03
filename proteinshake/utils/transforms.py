
class Compose():

    def __init__(self, *transforms):
        self.transforms = transforms

    def __call__(self, *args):
        for transform in self.transforms:
            if isinstance(args, (list, tuple)):
                data = transform(*args)
            else:
                data = transform(args)
        return data

    def __repr__(self):
        args = [f'  {transform}' for transform in self.transforms]
        return '{}([\n{}\n])'.format(self.__class__.__name__, ',\n'.join(args))
