from .containers import FunctionData
from ..population import Population
from operator import attrgetter

class Selection:
    """
    A class for handling all types of selection in the GA.
    
    Whenever a population needs to have some of its members selected
    for the creation of a parent pool or the next generation it 
    delegates this task to an instance of this class. The population has
    this instnace stored in its `ga_tools.selection` attribute.
    
    Each instance of this class supports being called. What a calling an
    instance does depends on the arguments the instance was initialized 
    with and what arguments were supplied during the call. In all cases
    the call implements a selection algorithm on a ``Population`` 
    instance.
    
    During initialization the functions which the ``Selection`` instance
    will use when its called are defined. Initialization therefore
    takes the names of methods defined in this class (as strings) and 
    saves them into the instance's `generational`, `mating` and 
    `mutation` attributes. These attributes should therefore always 
    hold the names of methods that are to be used for the given purpose
    - such as generational selection. The selection algorithms should be 
    written as methods within this class.
    
    When this class is called it requires a ``Population`` instance and 
    a string to be provided as arguments. The ``Population`` instance
    is the population which is to have some of its members selected.
    Consider the following code:
        
        >>> pop.select('generational')
    
    Here ``pop`` is a ``Population`` instance with a ``GATools`` 
    attribute, which holds an initialized ``Selection`` instance.
    
    The method `select` invoked in the code automatically provides
    the instance ``pop`` to its ``Selection`` instance when it carries 
    out a call to this instance. This means that each time `select` is 
    called on a population, the ``Selection`` instance will always act 
    on the population it is held by. Different populations can use
    different selection algorithms by holding ``Selection`` instances
    initialized with different method names. (Note that a ``Population``
    instance holds a ``Selection`` instance indirectly via its 
    `ga_tools` attribute.)
    
    The string provided to the `select` method is passed all the way
    down to the ``Selection`` instance being called. This is the second
    argument a ``Selection`` instance requires when it is being called.
    This is the string `generational` in the code example above. The
    string should be the name of one of the attributes of the 
    ``Selection`` class. This means that 'generational', 'mating' and
    'mutation' are valid strings at the time of this being written. If 
    more types of selection are added to MMEA, an attribute named after
    that type of selection should be added to the ``Selection`` class.
    If one wishes to invoke that type of selection on a population, 
    `select` must be called with the name of that attribute as a 
    parameter.
    
    There was a slight simplifcation in the paragraph describing 
    initialization. When the ``Selection`` instance is initialzed it
    is not just with the names of the selection methods to be used. It
    provided with the names of the methods and any paramters that the
    methods will need to use. For example the number of individuals
    that need to be selected. These parameters are packed into the
    ``FunctionData`` class. The ``FunctionData`` instance holding the
    method name and the appropriate parameter values is passed to
    the initializer of ``Selection``.
    
    Finally, because selection algorithms added as method to this class
    only need to be within the class but do not need access to the 
    instance of the class, the `self` parameter of methods is not 
    needed. This means the selection methods can be decorated with 
    ``staticmethod``.
    
    Attributes
    ----------
    generational : FunctionData

    mating : FunctionData

    mutation : FunctionData    
    
    """
    
    def __init__(self, generational, mating, mutation):
        """
        
        """
        
        self.generational = generational
        self.mating = mating
        self.mutation = mutation
    
    @classmethod
    def default(cls):
        func_data = FunctionData('fittest', size=5)
        return cls(*[func_data for x in range(0,3)])
    
    def __call__(self, population, type_):
        func_data = self.__dict__[type_]
        func = getattr(self, func_data.name)
        return func(population, **func_data.params)        

    @staticmethod
    def fittest(population, size):        
        if len(population) < size:
            raise ValueError(("Size of selected population" 
                              " must be less than or equal to" 
                              " size of original population."))                           
        
        ordered_pop = list(population.all_members())
        ordered_pop.sort(key=attrgetter('fitness'), reverse=True)    
        return Population(population.ga_tools, *ordered_pop[:size])
    
    @staticmethod    
    def roulette(population):
        pass
    @staticmethod
    def all_combinations(population):
        pass
    
    
    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """
