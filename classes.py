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
        
class Settings:
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
        return Population(*ordered_pop[:size], population.settings)
        
    def roulette(self, population):
        pass
    
    def all_combinations(self, population):
        pass

class Mating:
    def __init__(self, func_data):
        self.func_data = func_data
    
    def __call__(self, population):
        parent_pool = population.select('mating')
        offsprings = Population(population.settings)
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
    def __init__(self, *args):
        self.populations = []
        self.members = []
    
        for arg in args:
            if type(arg) is Population:
                self.populations.append(arg)
                continue
            
            if type(arg) is Cage:
                self.members.append(arg)
                continue
            
            if type(arg) is Settings and not hasattr(self, 'settings'):              
                self.settings = arg
                continue
                             
            raise TypeError(("Population must be"
                             " initialized with 'Population',"
                             " 'Cage' and 1 'Settings' type."))
                                    
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
    
    def select(self, type_):
        return self.settings.selection(self, type_)
        
    def gen_offspring(self):
        return self.settings.mating(self)
        
    def gen_mutants(self):
        return self.settings.mutation(self)
        
    def __iter__(self):
        return iter(self.all_members())
            
    def __getitem__(self, key):
        if type(key) is int:
            return list(self.all_members())[key]
        if type(key) is slice:
            cages = itertools.islice(self.all_members(), 
                                     key.start, key.stop, key.step)
            return Population(*cages, self.settings)
        
        raise ValueError("Index must be an integer or slice.")
        
    def __len__(self):
        return len(list(self.all_members()))
        
    def __sub__(self, other):        
        new_pop = Population(self.settings)        
        new_pop.add_members(*list(cage for cage in self 
                                                if cage not in other))
        return new_pop
        
    def __add__(self, other):
        return Population(self, other, self.settings)

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
                      Settings.default())
            pops.append(pop)
        
        pops[1].populations.append(pops[3])
        pops[1].populations.append(pops[4])
        pops[1].populations.append(pops[5])
        
        pops[2].populations.append(pops[6])
        pops[2].populations.append(pops[7])
        
        pops[0].populations.append(pops[1])
        pops[0].populations.append(pops[2])        
        
        return pops[0]
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        