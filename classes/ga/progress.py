import numpy as np

from ..population import Population


class GAProgress:
    def __init__(self):
        self.gens = []
        self.means = []
        self.mins = []
        self.maxs = []
        self.pops = Population()

    def update(self):
        self.gens.append(len(self.gens))
        self.past_pops.add_subpopulation(self.pop)
        
        if any(x.progress_params for x in self.pop):
            unscaled_var_mat = np.matrix([
                x.progress_params for x in self.pop if not 
                x.fitness_fail])

            self.maxs.append(np.max(unscaled_var_mat, axis=0).tolist()[0])
            self.mins.append(np.min(unscaled_var_mat, axis=0).tolist()[0])
            self.means.append(np.mean(unscaled_var_mat, axis=0).tolist()[0])
            
        else:
            self.means.append(self.pop.mean(lambda x : x.fitness))
            self.maxs.append(max(x.fitness for x in population))
            self.mins.append(min(x.fitness for x in population))
        

        