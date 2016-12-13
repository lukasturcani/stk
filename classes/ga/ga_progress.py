import numpy as np
      
class GAProgress:
    def __init__(self):
        self.gens = []
        self.means = []
        self.mins = []
        self.maxs = []

    def update(self, population):
        self.gens.append(len(self.gens))
        if hasattr(population[0], '_unscaled_fitness_vars'):
            unscaled_var_mat = np.matrix(
            [x.unscaled_fitness_vars.tolist() for 
             x in population if x.unscaled_fitness_vars is not None])

            self.maxs.append(np.max(unscaled_var_mat, axis=0).tolist()[0])
            self.mins.append(np.min(unscaled_var_mat, axis=0).tolist()[0])
            self.means.append(np.mean(unscaled_var_mat, axis=0).tolist()[0])
            
        else:
            self.means.append(population.mean(lambda x : x.fitness))
            self.maxs.append(max(x.fitness for x in population))
            self.mins.append(min(x.fitness for x in population))
        

        