import os
import numpy as np

class Mutation:
    def __init__(self, func_data):
        self.func_data = func_data
        self.n_calls = 0
        self.name = "mutation_{0}.mol"
    def __call__(self, population):
        self.n_calls += 1
        
        parent_pool = population.select('mutation')
        mutant_pop = Population(population.ga_tools)
        func = getattr(self, self.func_data.name)
        
        for parent in parent_pool:
            try:
                mutant = func(parent, **self.func_data.params)
                mutant_pop.members.append(mutant)
            except:
                continue

        mutant_pop -= population
            
        return mutant_pop


    def random_bb(self, cage, database):
        bb_file = np.random.choice(os.listdir(database))
        bb_file = os.path.join(database, bb_file)
        bb = BuildingBlock(bb_file)
        lk = next(x for x in cage.building_blocks if 
                                                isinstance(x, Linker))        
        return Cage((bb, lk), type(cage.topology), 
                                self.name.format(self.n_calls))
        
from ..population import Population
from ..molecular import BuildingBlock, Linker, Cage, Polymer