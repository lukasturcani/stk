import numpy as np
from functools import partial
from multiprocessing import Pool

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
        macro_mol.fitness = func(macro_mol)

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