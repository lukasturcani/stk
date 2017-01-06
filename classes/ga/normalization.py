from functools import partial
import numpy as np
import sys

class Normalization:
    def __init__(self, func_data):
        self.scaling_func = partial(getattr(self, func_data.name),
                                    **func_data.params)
                                    
    def __call__(self, population):
        self.scaling_func(population)        
        
    @staticmethod
    def carrots_and_sticks(population, carrot_coeffs, stick_coeffs,
                           carrot_exponents, stick_exponents):
        
        unscaled_carrots = [x.unscaled_fitness[0] for x in population if
                            isinstance(x.unscaled_fitness, tuple)]
                            
        unscaled_sticks = [x.unscaled_fitness[1] for x in population if
                           isinstance(x.unscaled_fitness, tuple)]
                           
        _carrot_means = np.mean(unscaled_carrots, axis=0)
        carrot_means = []
        for x in _carrot_means:
            if x == 0:
                carrot_means.append(1)
            else:
                carrot_means.append(x)
        
        _stick_means = np.mean(unscaled_sticks, axis=0)
        stick_means = []
        for x in _stick_means:
            if x == 0:
                stick_means.append(1)
            else:
                stick_means.append(x)
        
        for macro_mol in population:        
        
            # If one or more of the fitness parameters failed, 
            # return minimum fitness. 
            if macro_mol.fitness_fail:
                macro_mol.fitness = 1e-4
                continue
        
            # Calculate the scaled fitness parameters by dividing the 
            # unscaled ones by the fitness.
            scaled_carrots = np.divide(macro_mol.unscaled_fitness[0], 
                                       carrot_means)
              
            scaled_sticks = np.divide(macro_mol.unscaled_fitness[1],
                                      stick_means)        
              
            try:  
                scaled_carrots = np.power(scaled_carrots, 
                                          carrot_exponents)
                scaled_sticks = np.power(scaled_sticks, 
                                         stick_exponents)
                                         
                scaled_carrots = np.multiply(scaled_carrots, 
                                             carrot_coeffs)
                scaled_sticks = np.multiply(scaled_sticks, 
                                            stick_coeffs)

            # If the user forgot put the wrong number values in the 
            # exponent or coefficient arrays MMEA will tell them and 
            # exit gracefully.
            except  ValueError:
                print(('Fitness function calculates carrot array of '
                       'size {} and stick array of size {}. This does '
                       'not match the array sizes given to the '
                       'normalization function:\n'
                       '\tcarrot_coeffs : {}\n'
                       '\tcarrot_exponents : {}\n'
                       '\tstick_coeffs : {}\n'
                       '\tstick_exponents : {}\n').format(
                                   len(scaled_carrots),
                                   len(scaled_sticks),
                                   len(carrot_coeffs),
                                   len(carrot_exponents),
                                   len(stick_coeffs),
                                   len(stick_exponents)))
                sys.exit()

            carrot_term = np.sum(scaled_carrots)
            penalty_term = np.sum(scaled_sticks)
            penalty_term =  np.divide(1,penalty_term)
            if penalty_term > 1e101:
                penalty_term = 1e101
                        
            macro_mol.fitness = penalty_term + carrot_term    