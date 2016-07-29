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
