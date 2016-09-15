import numpy as np
from functools import partial
from multiprocessing import Pool
import itertools
import rdkit.Chem as chem
from scipy.spatial.distance import euclidean

from .pyWindow import window_sizes
from .classes.exception import MacroMolError

def calc_fitness(func_data, population):
    """
    Calculates the fitness values of all members of a population.
    
    A fitness function should take a ``MacroMolecule`` instance and
    return a number representing its fitness. The assignement to the
    `fitness` attribute of a population member happens here, not by the
    fitness function.    
    
    Parameters
    ----------
    func_data : FunctionData
        A ``FunctionData`` instance representing the chosen fitness 
        function and any additional parameters it may require.
    
    population : Population
        The population whose members must have their fitness calculated.
        
    Returns
    -------
    None : NoneType    
    
    """

    # Get the fitness function object.
    func = globals()[func_data.name]

    # Apply the function to every member of the population.
    for macro_mol in population:
        try: 
            macro_mol.fitness = func(macro_mol, **func_data.params)
            
        except Exception as ex:
            MacroMolError(ex, macro_mol, 'During fitness calculation.')
            macro_mol.fitness = -np.inf
            macro_mol.windows = None

        print(macro_mol.fitness, '-', macro_mol.prist_mol_file)            
            
def random_fitness(macro_mol):
    """
    Returns a random fitness value between 0 and 100.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to which a fitness value is to be assigned.
    
    Returns
    -------
    int
        An integer between 0 (including) and 100 (excluding).

    """

    return np.random.randint(0,100)

def cage(macro_mol, target_size, coeffs=None, exponents=None):
    if macro_mol.fitness:
        print('Skipping {0}'.format(macro_mol.prist_mol_file))
        return macro_mol.fitness
    
    if macro_mol.topology.windows is None:
        return -np.inf

    if coeffs is None:
        coeffs = np.array([50,1,1])
        
    if exponents is None:
        exponents = np.array([1,1,1])

    if target_size <= 10:
        coeffs[0] = 5
        exponents[0] = 2        

    cavity_diff = abs(target_size - macro_mol.topology.cavity_size())

    target_window_area = np.square(target_size)
    window_area = np.square(max(macro_mol.topology.windows))
    window_area_diff = abs(target_window_area - window_area)
            
    fitness_value = np.array([
                             cavity_diff, 
                             window_area_diff,                                                          
                             macro_mol.topology.window_difference(500)
                             ])
    
    fitness_value = np.power(fitness_value, exponents)
    fitness_value = np.multiply(fitness_value, coeffs)    

    return -np.sum(fitness_value)

    
    