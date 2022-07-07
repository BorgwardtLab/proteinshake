import os
import torch

# decorator to return saved file if it exists
def checkpoint(path):
    def decorator(function):
        def wrapper(self, *args, **kwargs):
            if os.path.exists(path.format(**self.__dict__)):
                return torch.load(path.format(**self.__dict__))
            else:
                os.makedirs(os.path.dirname(path.format(**self.__dict__)), exist_ok=True)
                result = function(self, *args, **kwargs)
                torch.save(result, path.format(**self.__dict__))
                return result
        return wrapper
    return decorator
