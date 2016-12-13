from functools import partial
import numpy as np

class Normalization:
    def __init__(self, func_data):
        self.scaling_func = partial(getattr(self, func_data.name),
                                    **func_data.params)
                                    
    def __call__(self, population):
        self.scaling_func(population)        
        
    @staticmethod
    def cage(population):
        # If one or more of the fitness parameters failed, return minimum 
        # fitness. 
        if macro_mol._fitness_fail:
            return 1e-4
    
        # Set the default coeffient values.
        if coeffs is None:
            coeffs = np.array([1,1,1,1,0.2])
            
        # Set the default exponent values.
        if exponents is None:
            exponents = np.array([1,1,1,1,1]) 
    
        # Calculate the scaled fitness parameters by dividing the unscaled
        # ones by the fitness.
        scaled = np.divide(macro_mol._unscaled_fitness_vars, 
                           population._fitness_cache['mean'])
           
        fitness_vars = np.power(scaled, exponents)
        fitness_vars = np.multiply(fitness_vars, coeffs)    
        penalty_term = np.sum(fitness_vars[:-1])
        penalty_term =  np.divide(1,penalty_term)
        if penalty_term > 1e101:
            penalty_term = 1e101
        
        # Carrots and sticks, where the previous fitness parameters were
        # the sticks.
        carrot_term = fitness_vars[-1]
        
        return penalty_term + carrot_term    