import numpy as np
from functools import wraps
from operator import attrgetter
import itertools

class Cage:
    def __init__(self):
        pass  

    def same_cage(self, other):
        return self.bb == other.bb and self.lk == other.lk
        
    def __str__(self):
        return str(self.__dict__) + "\n"
    
    def __repr__(self):
        return str(self.__dict__) + "\n"
    
    

    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """

    @classmethod
    def init_empty(cls):
        obj = cls()
        string = ['a','b','c','d','e','f','g','h','i','j','k','l','m',
                  'n','o', 'p','q','r','s','t','u','v','w','x','y','z']
        obj.bb = np.random.choice(string)
        obj.lk = np.random.choice(string)
        obj.fitness = abs(np.random.sample())
        return obj

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

class Selection:
    def __init__(self, generational, mating, mutation):
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

    def fittest(self, population, size):        
        if len(population) < size:
            raise ValueError(("Size of selected population" 
                              " must be less than or equal to" 
                              " size of original population."))                           
        
        ordered_pop = list(population.all_members())
        ordered_pop.sort(key=attrgetter('fitness'), reverse=True)    
        return Population(*ordered_pop[:size], population.ga_tools)
        
    def roulette(self, population):
        pass
    
    def all_combinations(self, population):
        pass

class Mating:
    def __init__(self, func_data):
        self.func_data = func_data
    
    def __call__(self, population):
        parent_pool = population.select('mating')
        offsprings = Population(population.ga_tools)
        func = getattr(self, self.func_data.name)
        
        for parents in parent_pool.populations:
            offspring = func(*parents, **self.func_data.params)
            offsprings.add_members(offspring)

        offsprings -= population
            
        return offsprings
        
    def bb_lk_exchange(self, cage1, cage2, _):
        ...

class Mutation:
    def __init__(self):
        pass

class Population:
    """
    A container for instances of ``Cage`` and ``Population``.

    This is the central class of MMEA. The GA is invoked by calling the
    ``gen_offspring``, ``gen_mutants`` and ``select`` methods of this 
    class on a given instance. However, this class is a container for ``Cage`` and other 
    ``Population`` instances first and foremost. It delegates GA operations to
    its `ga_tools` attribute. Any functionality related to the GA should be 
    delegated to this attribute, as done with the aforementioned methods. 
    A comphrehensive account of how this delegation is implemented is provided
    in the module level docstring.
    
    For consistency and maintainability, collections ``Cage`` or ``Population`` instances should always be placed
    in a ``Population`` instance. As a result, any function which should return
    multiple ``Cage`` or ``Population`` instances can be expected
    to return a single ``Population`` instance holding the desired instances.     
    
    The only operations directly addressed by this class and definined within
    it are those relevant to its role as a container. It supports
    all expected and necessary container operations such as iteration, indexing, membership 
    checks (via the 'is in' operator) as would be expected. Additional 
    operations such as comparison via the ``==``, ``>``, etc. operators is also supported. 
    Details of the various implementations and a full list of supported operations 
    can be found by examining the included methods. Note that all comparison
    operations are accounted for with the ``total_ordering`` decorator, even
    if they are not explicity defined.

    Attributes
    ----------
    populations : list
        A list of other instances of the ``Population`` class. This
        allows the implementation of subpopulations or 'evolutionary 
        islands'. This attribute is also used for grouping cages within
        a given population such as when grouping parents together from 
        within a parent pool.
        
    members : list
        A list of ``Cage`` instances. These are the members of the
        population which are not held within any subpopulations. This 
        means that not all members of a population are stored here. 
        To access all members of a population the generator method 
        ``all_members`` should be used.
    
    ga_tools : GATools, optional
        An instance of the ``GATools`` class. Calls to preform GA
        operations on the ``Population`` instance are delegated 
        to this attribute.
    
    """
    
    def __init__(self, *args):
        """
        Initializer for ``Population`` class.
        
        This initializer creates a new population from the 
        ``Population`` and ``Cage`` instances given as arguments.
        It also accepts a single, optional ``GATools`` instance if the 
        population is to have GA operations performed on it.
        
        The arguments can be provided in any order regardless of type.
        
        Parameters
        ----------
        *args : Cage, Population, GATools
            A population is initialized with as many ``Cage`` or 
            ``Population`` arguments as required. These are placed into 
            the `members` or `populations` attributes, respectively. 
            A single ``GATools`` instance may be included and will be
            placed into the `ga_tools` attribute. 
        
        Raises
        ------
        TypeError
            If the instance is initialized with something other than
            ``Cage``, ``Population`` or more than 1 ``GATools`` object. 
            
        """    
        
        self.populations = []
        self.members = []
    
        for arg in args:
            if type(arg) is Population:
                self.populations.append(arg)
                continue
            
            if type(arg) is Cage:
                self.members.append(arg)
                continue
            
            if type(arg) is GATools and not hasattr(self, 'ga_tools'):              
                self.ga_tools = arg
                continue
                             
            raise TypeError(("Population must be"
                             " initialized with 'Population',"
                             " 'Cage' and 1 'GATools' type."))
                                    
    def all_members(self):       
        for ind in self.members:
            yield ind
        
        for pop in self.populations:
            yield from pop.all_members()
            
    def add_members(self, *args, duplicates=False):
        for arg in args:
            if type(arg) is Population:
                if not duplicates:
                    self.members.append(cage for cage in arg 
                                                    if cage not in self)
                if duplicates:
                    self.members.append(cage for cage in arg)

            if type(arg) is Cage:
                is_duplicate = arg in self
                if duplicates:
                    self.members.append(arg)
                if not duplicates and not is_duplicate:
                    self.members.append(arg)                                    
    
    def select(self, type_='generational'):
        return self.ga_tools.selection(self, type_)
        
    def gen_offspring(self):
        return self.ga_tools.mating(self)
        
    def gen_mutants(self):
        return self.ga_tools.mutation(self)
        
    def __iter__(self):
        return iter(self.all_members())
            
    def __getitem__(self, key):
        if type(key) is int:
            return list(self.all_members())[key]
        if type(key) is slice:
            cages = itertools.islice(self.all_members(), 
                                     key.start, key.stop, key.step)
            return Population(*cages, self.ga_tools)
        
        raise ValueError("Index must be an integer or slice.")
        
    def __len__(self):
        return len(list(self.all_members()))
        
    def __sub__(self, other):        
        new_pop = Population(self.ga_tools)        
        new_pop.add_members(*list(cage for cage in self 
                                                if cage not in other))
        return new_pop
        
    def __add__(self, other):
        return Population(self, other, self.ga_tools)

    def __contains__(self, item):
        return any(item.same_cage(cage) for cage in self.all_members())

    def __str__(self):
        output_string = (" Population " + str(id(self)) + "\n" + 
                            "--------------------------\n" + 
                            "\tMembers\n" + "   ---------\n")
        
        for cage in self.members:
            output_string += "\t"  + str(cage) + "\n"
        
        output_string += (("\tSub-populations\n" + 
                           "   -----------------\n\t"))
        
        for pop in self.populations:
            output_string += str(id(pop)) + ", "
        
        output_string += "\n\n"

        for pop in self.populations:
            output_string += str(pop)
        
        return output_string
        
    def __repr__(self):
        return str(self)
        


    """
    The following methods are inteded for convenience while 
    debugging or testing and should not be used during typical 
    execution of the program.
    
    """

    @classmethod
    def init_empty(cls):
        pops = []
        
        for x in range(0,8):
            pop = cls(*iter(Cage.init_empty() for x in range(0,3)), 
                      GATools.default())
            pops.append(pop)
        
        pops[1].populations.append(pops[3])
        pops[1].populations.append(pops[4])
        pops[1].populations.append(pops[5])
        
        pops[2].populations.append(pops[6])
        pops[2].populations.append(pops[7])
        
        pops[0].populations.append(pops[1])
        pops[0].populations.append(pops[2])        
        
        return pops[0]
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        