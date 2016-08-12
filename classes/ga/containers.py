class FunctionData:
    """
    Stores information about functions and their parameters.

    Attributes
    ----------
    name : str
        The name of a function or method.

    params : dict
        The parameters of the function or method who's name is held by
        'name' which should be used when that function is called.
    
    """
    
    __slots__ = ['name', 'params']    
    
    def __init__(self, name, **kwargs):
        self.name = name
        self.params = kwargs
        
class GATools:
    """
    Stores ``Selection``, ``Mating`` and ``Mutation`` instances.
    
    Instances of this class are held by ``Population`` instances in
    their `ga_tools` attribute. All the instances which carry out GA
    operations on the population are held conveniently together in an 
    instance of this class.
    
    This class also holds the name of the method in the ``Optimization``
    class which is to be used for optimizing the structures of
    ``MacroMolecule`` instances.

    Attributes
    ----------
    selection: Selection
        The ``Selection`` instance which performes selections of a 
        population's members.
    
    mating : Mating
        The ``Mating`` instance which mates a population's members.
    
    mutation : Mutation
        The ``Mutation`` instance which mutates a population's members
        
    opt_func_name : str
        The name of the method in the ``Optimization`` class to be used
        for ``MacroMolecule`` optimizations.
    
    """
    
    __slots__ = ['selection', 'mating', 'mutation', 'optimization']    
    
    def __init__(self, selection, mating, mutation, optimization):
        self.selection = selection
        self.mating = mating
        self.mutation = mutation
        self.optimization = optimization

