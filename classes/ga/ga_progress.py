import numpy as np
import os
import matplotlib.pyplot as plt
        
class GAProgress:
    def __init__(self):
        self.gens = []
        self.means = []
        self.mins = []
        self.maxs = []

    def update(self, population):
        self.gens.append(len(self.gens))
        if hasattr(population[0], 'unscaled_fitness_vars'):
            unscaled_var_mat = np.matrix(
            [x.unscaled_fitness_vars.tolist() for 
             x in population if x.unscaled_fitness_vars is not None])
            
            self.maxs.append(np.max(unscaled_var_mat, axis=0).tolist()[0])
            self.mins.append(np.min(unscaled_var_mat, axis=0).tolist()[0])
            print(self.means)
            
        else:
            self.means.append(population.mean('fitness'))
            self.maxs.append(max(x.fitness for x in population))
            self.mins.append(min(x.fitness for x in population))
        
    def epp(self, plot_name):
        if (isinstance(self.maxs[0], float) or 
                                        isinstance(self.maxs[0], int)):
            fig = plt.figure()
            plt.plot(self.gens, self.means, color='green')
            plt.plot(self.gens, self.mins, color='blue')
            plt.plot(self.gens, self.maxs, color='red')
            fig.savefig(plot_name, dpi=fig.dpi)
            plt.close('all')
            
        else:
            for x in range(len(self.means[0])):
                fig = plt.figure()
                y_mean = [v[x] for v in self.means]
                y_max = [v[x] for v in self.maxs]
                y_min = [v[x] for v in self.mins]

                plt.plot(self.gens, y_mean, color='green')
                plt.plot(self.gens, y_min, color='blue')
                plt.plot(self.gens, y_max, color='red')
                
                new_plot_name = str(x).join(os.path.splitext(plot_name))                
                
                fig.savefig(new_plot_name, dpi=fig.dpi)
                plt.close('all')

            
        