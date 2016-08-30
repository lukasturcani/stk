import numpy as np
from functools import partial
from multiprocessing import Pool

def calc_fitness(func_data, population):
    func = globals()[func_data.name]

    # Apply the function to every member of the population.
    
    for macro_mol in population:
        macro_mol.fitness = func(macro_mol)
    

def random_fitness(macro_mol):
    return np.random.randint(0,100)