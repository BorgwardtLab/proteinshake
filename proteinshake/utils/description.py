
class DatasetDescription:

    def __init__(self, annotation):
        self.Annotation = annotation

    def to_dict(self):
        return {key:getattr(self, key) for key in []}
    

class TaskDescription:

    def __init__(self, type, input, output):
        self.Type = type
        self.Input = input
        self.Output = output

    def to_dict(self):
        return {key:getattr(self, key) for key in []}