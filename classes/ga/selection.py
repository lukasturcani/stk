from operator import attrgetter
import itertools
 
# Note that import of the Population class occurs at the end of the
# file. There is a good reason for this.
from .containers import FunctionData

class Selection:
    """
    A class for handling all types of selection in the GA.
    
    Whenever a population needs to have some of its members selected
    for the creation of a parent pool or the next generation it 
    delegates this task to an instance of this class. The population has
    this instance stored in its `ga_tools.selection` attribute.
    
    Each instance of this class supports being called. What a calling an
    instance does depends on the arguments the instance was initialized 
    with and what arguments were supplied during the call. In all cases
    the call implements returns a generator. This generator yields
    members of a ``Population`` instance in accordance to some selection
    algorithm.
    
    Initialization of this class takes the names of methods defined in 
    this class (as strings) and saves them into the instance's 
    `generational`, `mating` and `mutation` attributes. These attributes 
    should therefore always hold the names of methods that are to be 
    used for the given purpose - such as generational selection. The 
    selection algorithms should be written as methods within this class.
    
    When an instance of this class is called it requires a 
    ``Population`` instance and a string to be provided as arguments. 
    The ``Population`` instance is the population which is to have some 
    of its members selected. Consider the following code:
        
        >>> pop.select('generational')
    
    Here ``pop`` is a ``Population`` instance with a `ga_tools` 
    attribute, which holds an initialized ``Selection`` instance.
    
    The method `select` invoked in the code automatically provides
    the instance ``pop`` to the ``Selection`` instance it contains. The 
    function then calls the ``Selection`` instance. This means that each 
    time `select` is called on a population, the ``Selection`` instance 
    will always act on the population it is held by. Different 
    populations can use different selection algorithms by holding 
    different ``Selection`` instances. Alternatively, they can perform
    the same selection algorithms by sharing a ``Selection`` instance.
    
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
    methods will need to use. These parameters are packed into the 
    ``FunctionData`` class. The ``FunctionData`` instance holding the
    method name and the appropriate parameter values is passed to
    the initializer of ``Selection``.
    
    Selection algorithms added as methods to this class only need to be 
    within the class but do not need access to the instance of the 
    class, the `self` parameter of methods is not needed. This means the 
    selection methods can be decorated with ``staticmethod``. In the 
    cases where access to the class or instance is necessary, 
    ``classmethod`` or no decorator should be used as appropriate.
    
    Selection algorithms are to be implemented as generators. Selection 
    algorithms which produce parents pools must yield a tuple of
    ``MacroMolecule`` instances. Selection algorithms should be grouped
    together by their expected use when written into the class body.    
    
    Attributes
    ----------
    generational : FunctionData
        This holds the ``FunctionData`` object representing the
        selection function used for selecting the next generation of 
        individuals along with any parameters the function may require.

    mating : FunctionData
        Holds the ``FunctionData`` object representing the selection 
        function used for selecting parents, along with any parameters
        the function may require.
        
    mutation : FunctionData    
        Holds the ``FunctionData`` object representing the selection
        function used for selection individuals for mutation, along with
        any parameters the function may require.
    
    """
    
    def __init__(self, generational, mating, mutation):        
        self.generational = generational
        self.mating = mating
        self.mutation = mutation
    
    def __call__(self, population, type_):
        """
        Implements the selection algorithm chosen during initializaiton.        
        
        Parameters
        ----------
        population : Population
            The population from which members should be selected.

        type_ : str
            The name of a ``Selection`` attribute. The name corresponds
            to the type of selection that is desired, ie generational,
            mating or mutation.
            
        Returns
        -------
        generator
            A generator which yields selected members of the population.
            The generator will yield ``MacroMolecule`` instances unless
            it yields parents for mating. In this case it yields a tuple
            of ``MacroMolecule`` instances.
        
        """
        
        # Get the attribute with the name provided as a string in 
        # `type_`. This returns a ``FunctionData`` object holding the
        # name of the method which needs to be implemented and any
        # required parameters and their values.
        func_data = self.__dict__[type_]
        # Using the name of the method, get the method object with 
        # ``getattr``.
        func = getattr(self, func_data.name)
        # Call the method on the population and provide any additional
        # parameters which may be necessary.
        return func(population, **func_data.params)        

    """
    The following selection algorithms can be used for generational
    selection and the selection of mutants. They cannot be used for
    selection of parents.
    
    """

    @staticmethod
    def fittest(population):        
        """
        Yields members of the population, fittest first.
        
        This yields members regardless of which subpopulation they are
        in.        
        
        Parameters
        ----------
        population : Population
            The population who's members need to be selected.
        
        Yields
        ------
        MacroMolecule
            The next fittest ``MacroMolecule`` instance held by the
            population.
        
        """        
          
        for ind in sorted(Population, reverse=True):
            yield ind
    
    @staticmethod    
    def roulette(population):
        pass

    """
    The following selection algorithms can be used for the selection of
    parents only.    
    
    """

    @staticmethod
    def all_combinations(population):
        """
        Yields every possible pairing of individuals from a population.
        
        This yields members regardless of which subpopulation they are 
        in. Each pair is only returned once. This means if (1,2) is 
        returned (2,1) will not be.
        
        Parameters
        ----------
        population : Population
            The population from which parents should be selected.        
        
        Yields
        ------
        tuple of 2 MacroMolecule instances
            The ``MacroMolecule`` instances which together form a parent
            pair.
        
        """
        
        # Get all combinations of size 2.        
        for mol1, mol2 in itertools.combinations(population, 2):
            yield mol1, mol2 
    
    @classmethod
    def all_combinations_n_fittest(cls, population, n):
        """
        Yields all pairings of the `n` fittest individuals.

        Parameters
        ----------
        population : Population
            The population from which parents should be selected.
            
        Yields
        ------
        tuple of 2 MacroMolecule instances
            The ``MacroMolecule`` instances which together form a parent
            pair.
        
        """
        
        n_fittest = itertools.islice(cls.fittest(population), n)
        for ind1, ind2 in itertools.combinations(n_fittest, 2):
            yield ind1, ind2
    
    
    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """

    @classmethod
    def default(cls):
        func_data = FunctionData('fittest', size=5)
        return cls(*[func_data for x in range(0,3)])
        
        
from ..population import Population
