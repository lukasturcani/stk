"""
Defines the FunctionData class.

"""

class FunctionData:
    """
    Stores information about functions and their parameters.

    Attributes
    ----------
    name : str
        The name of a function or method.

    params : dict
        The parameters of the function or method who's name is held by
        `name`.
    
    """
    
    __slots__ = ['name', 'params']    
    
    def __init__(self, name, **kwargs):
        self.name = name
        self.params = kwargs
        
    def __hash__(self):
        return 1

    def __eq__(self, other):
        
        same_length = len(self) == len(other)
        same_items = all(x in other.params.items() for x in self.params.items())        
        same_name = self.name == other.name        
        
        return same_length and same_items and same_name
        
    def __len__(self):
        return len(self.params.items())        
        
    def __str__(self):
        s = ", ".join("{}={!r}".format(key, value) for key, value in self.params.items())
        return "FunctionData({!r}, ".format(self.name) + s + ")"
    def __repr__(self):
        return str(self)