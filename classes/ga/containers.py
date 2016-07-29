class FunctionData:
    def __init__(self, name, **kwargs):
        self.name = name
        self.params = kwargs
        
class GATools:
    def __init__(self, selection, mating, mutation):
        self.selection = selection
        self.mating = mating
        self.mutation = mutation
        
    @classmethod
    def default(cls):
        return cls(Selection.default(),2,3)