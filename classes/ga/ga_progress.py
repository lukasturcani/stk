import numpy as np
      
class GAProgress:
    def __init__(self):
        self.gens = []
        self.means = []
        self.mins = []
        self.maxs = []

    def update(self, population):
        self.gens.append(len(self.gens))
        if any(x.progress_params for x in population):
            unscaled_var_mat = np.matrix([
                x.progress_params for x in population if not 
                x.fitness_fail])

            self.maxs.append(np.max(unscaled_var_mat, axis=0).tolist()[0])
            self.mins.append(np.min(unscaled_var_mat, axis=0).tolist()[0])
            self.means.append(np.mean(unscaled_var_mat, axis=0).tolist()[0])
            
        else:
            self.means.append(population.mean(lambda x : x.fitness))
            self.maxs.append(max(x.fitness for x in population))
            self.mins.append(min(x.fitness for x in population))
        

        