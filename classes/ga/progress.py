"""
Defines the GAProgress class.

"""

import numpy as np
from operator import attrgetter
import copy

from ..population import Population


class GAProgress:
    """
    Stores the GA's data across generations.

    This class is used for tracking how the fitness values and other
    parameters of population members change across generations.
    
    Some fitness functions may calculate values besides the fitness,
    which are desired for tracking, for example energy. If this is the
    case, they place the energy in the `progress_params` attribute of 
    the MacroMolecule. These `progress_params` are then tracked here.
    
    For more information on `progress_params` see the documentation in
    ``mmea/fitness.py``.
    
    Attributes
    ----------
    gens : list of ints
        The numbers of generations.
    
    means : list of floats or of lists of floats ie 
            [[float, float, ...], [float, float, ..], ... ]
            
        This list holds the mean fitness at each genertion in `gens`.
        Alternatively, if `progress_params` was used by the fitness
        function, the list holds the mean values of each element in 
        `progress_params` at each generation.
    
    mins : list of floats or of lists of floats ie 
           [[float, float, ...], [float, float, ..], ... ] 
           
        This list holds the minimum fitness at each genertion in `gens`.
        Alternatively, if `progress_params` was used by the fitness
        function, the list holds the minimum values of each element in 
        `progress_params` at each generation.            
    
    maxs : list of floats or of lists of floats ie 
           [[float, float, ...], [float, float, ..], ... ]
           
        This list holds the maximum fitness at each genertion in `gens`.
        Alternatively, if `progress_params` was used by the fitness
        function, the list holds the maximum values of each element in 
        `progress_params` at each generation.
        
    pops : Population
        This population holds each generations Population instance as 
        one of its subpopulations.
        
    hist : list of GAProgress instances
        Using the ``normalization()`` method overwrites the values of
        the `means`, `mins` and `maxs` attributes.  In order to keep
        the old values, each time this method is run, a copy of the 
        GAProgress instances with the old values is saved in this
        attribute.
        
    """
    
    def __init__(self):
        self.gens = []
        self.means = []
        self.mins = []
        self.maxs = []
        self.pops = Population()
        self.hist = []

    def update(self, pop):
        """
        Updates the attributes with values from `pop`.
        
        Parameters
        ----------
        pop : Population
            The population whose values are used to update the values 
            in the attributes of GAProgress.
        
        Returns
        -------
        None : NoneType
        
        """
        
        self.gens.append(len(self.gens))
        self.pops.add_subpopulation(pop)
        
        # If the MacroMolecules have the `progress_params` paramter
        # then track these. Otherwise track fitness values.
        if any(x.progress_params for x in pop):
            unscaled_var_mat = np.matrix([
                x.progress_params for x in pop if not 
                x.fitness_fail])

            self.maxs.append(np.max(unscaled_var_mat, 
                                    axis=0).tolist()[0])
            
            self.mins.append(np.min(unscaled_var_mat, 
                                    axis=0).tolist()[0])
            
            self.means.append(np.mean(unscaled_var_mat, 
                                      axis=0).tolist()[0])
            
        else:
            self.means.append(pop.mean(lambda x : x.fitness))
            self.maxs.append(max(x.fitness for x in pop))
            self.mins.append(min(x.fitness for x in pop))
        

    def normalize(self, norm_func):
        """
        Normalizes the fitness of members in `pops` and updates attrs.
        
        Applies the `norm_func` normalization function to the Population
        `pops`. It then uses the new fitness values to update the
        values in `means`, `mins` and `maxs`.
        
        
        Parameters
        ----------
        norm_func : functools.partial
            A normalization function.
        
        Modifies
        --------
        hist : list of GAProgress instances
            A copy of the GAProgress instance before the normalization
            is placed here.
        
        means : kist of floats
            The means now reflect the fitness values of each
            subpopulation in `pops`.
        
        mins : list of floats
            The minimums now reflect the fitness values of each
            subpopulation in `pops`.

        maxs : list of floats
            The maximums now reflect the fitness values of each
            subpopulation in `pops`.
        
        Returns
        -------
        None : NoneType
        
        """
        
        self.hist.append(copy.deepcopy(self))        
        
        norm_func(self.pops)
        
        self.gens = []
        self.mins = []
        self.means = []
        self.maxs = []
        for i, gen in enumerate(self.pops.populations, 1):
            self.gens.append(i)
            self.mins.append(min(gen, key=attrgetter('fitness')).fitness)
            self.means.append(gen.mean(attrgetter('fitness')))
            self.maxs.append(max(gen, key=attrgetter('fitness')).fitness)