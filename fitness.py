"""
Module for defining fitness functions.

Extending MMEA: Adding fitness functions
----------------------------------------
To add a new fitness function simply write it as a function in this
module. It will need to take the ``MacroMolecule`` instance as its first
argument and this argument should be called ``macro_mol``. The purpose
of this is to help users identify which arguments are handled
automatically by MMEA and which they need to define in the input file.
The convention is that if the fitness function takes an argument called
``macro_mol`` they do not have to specify that argument in the input 
file. 

Optionally, a fitness function may take an argument ``population``. If
this is done, it should be the argument after ``macro_mol``. When this
argument is present the population is also passed to the fitness
function, along with the MacroMolecule instance. This can be useful if
some sort of normalization is desired.

Fitness functions should return a value. This value represents the 
fitness of the MacroMolecule instance in the argument `macro_mol`. They 
do not assign the value to the ``MacroMolecule`` instance. This is done
by the ``_calc_fitness()`` function. Don't worry about it.

A fitness function may be complex and may not fit neatly into a single 
function. For example, the ``cage_target()`` fitness function needs to 
call ``_generate_complexes()`` in order to sample various conformations
before outputting a fitness value. This is fine. Define helper functions 
such as ``_generate_complexes()`` within this module but make sure they 
are private. This means that names of helper functions begin with a 
leading underscore. 

"""

import numpy as np
import rdkit.Chem as chem
import copy
from inspect import signature
from functools import partial

from .classes.exception import MacroMolError
from .classes.molecular import MacroMolecule, StructUnit
from . import optimization
from .convenience_tools import rotation_matrix_arbitrary_axis

def _calc_fitness(func_data, population):
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
    
    # Check if the fitness function wants to use the population as well.
    use_pop = 'population' in signature(func).parameters.keys()
    # Fitness functions can cache values here each generation. For
    # example, storing the mean values of cavity differences each
    # generation is used by the ``cage()`` fitness function.
    population._fitness_cache = {}

    # Apply the function to every member of the population.
    for macro_mol in population:
        try:
            if use_pop:
                macro_mol.fitness = func(macro_mol, population, 
                                         **func_data.params)
            else:
                macro_mol.fitness = func(macro_mol, **func_data.params)
                
        except Exception as ex:
            MacroMolError(ex, macro_mol, 'During fitness calculation.')
            raise ex

    # After each macro_mol has a fitness value, sort the population by 
    # fitness and print.
    for macro_mol in sorted(population, reverse=True):
        print(macro_mol.fitness, '-', macro_mol.prist_mol_file)

def _param_labels(*labels):
    """
    Adds `param_labels` attribute to a fitness function.
    
    Fitness functions which undergo the scaling procedure have an EPP 
    graph plotted for each attribute used to calculate total fitness. 
    For example, if the ``cage`` fitness function was used
    during the GA run 5 graphs would be plotted at the end of the run.
    One for each unscaled ``var`` (see ``cage`` documentation for more 
    details). This plotting is done by the ``GAProgress`` class. 
    
    In order for the ``GAProgress`` class to produce decent graphs, 
    which means that each graph has the y-axis and title labeled with 
    the name of the ``var`` plotted, it needs to know what each ``var`` 
    should be called.
    
    This is done by adding a `param_labels` attribute to the fitness
    function. The ``GAProgress`` class acccesses this attribute during
    plotting and uses it to set the y-axis / title.
    
    Parameters
    ----------
    labels : tuple
        List of strings about the fitness labels used for plotting EPPs.
        The order of the strings should represent the order of the
        fitness ``vars`` in the fitness funciton. In practice it should
        correspond to the order of the ``coeffs`` or ``exponents`` 
        parameters given to the fitness function.
    
    Returns
    -------
    func
        Decorated function.
        
    """
    
    def add_labels(func):
        func.param_labels = labels    
        return func
        
    return add_labels
            
def random_fitness(macro_mol):
    """
    Returns a random fitness value between 1 and 10.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to which a fitness value is to be assigned.
    
    Returns
    -------
    int
        An integer between 0 (including) and 100 (excluding).

    """

    return np.random.randint(1,10)

# Calls the decorator with the specific labels
@_param_labels('Cavity Difference ','Window Difference ',
                'Asymmetry ', 'Negative Energy per Bond ', 
                'Positive Energy per Bond ')
def cage(macro_mol, population, target_cavity,  
         target_window=None, coeffs=None, exponents=None, 
         energy_params={'key':('rdkit', 'uff')}):
    """
    Calculates the fitness of a cage.
    
    The fitness function has the form
    
        (1) fitness = penalty_term + carrot_term
        
    where
        
        (2) penalty_term = 1 / [ A*(var1^a) + B*(var2^b) + C*(var3^c) + 
                                 D*(var4^d) ],
                                 
    and 
    
        (3) carrot_term = E*(var5^e).
        
    Here var1 to var5 signify parameters of the cage which factor into
    fitness. These parameters are calculated by the fitness function and
    placed in variables, for example:
        
        1) `cavity_diff` - the difference between the cage's cavity
           (diameter) and the desired cavity size
        2) `window_diff` - the difference between the diameter of the 
           window of the cage and the size required to allow the target
           to enter the cavity
        3) `asymmetry` - sum of the difference between the size of the
           windows of the same type
        4) `neg_eng_per_bond` - the formation energy of the cage divided 
           by the number of bonds for during assembly, when the energy
           is < 0. This is a measure of the stability or strain of a 
           cage.
        5) Same as 4) but used when the energy is > 0.
        
    The design of the fitness function is as follows. Consider two
    cages, ``CageA`` and ``CageB``. If the parameters, 1) to 4), in
    CageA, which signify poor fitness are 10 times that of CageB and the
    parameter which signifies good fitness, 5), is 10 times less than 
    CageB, then CageA should have a fitness value 10 times less than 
    CageB.
    
    This is assuming all coefficients and powers are 1. Note that a cage
    will always have either 4) or 5) at 0. This is because its total 
    energy will always be either positive or negative.
    
    The `coeffs` parameter has the form
    
        np.array([1,2,3,4,5]),

    where 1, 2, 3, 4 and 5 correspond to the desired values of A, B, C,  
    D and E in equation (2). Equally the `exponents` parameter also has 
    the form

        np.array([5,6,7,8,9]),
    
    where 5, 6, 7, 8 and 9 correspond to the values of a, b, c, d and e
    in equation (2).
    
    Assume that for a given GA run, it is not worthwhile factoring in 
    the `window_diff` parameter. This may be because cages which 
    form via templating are also to be considered. If the cage forms
    around the target, its windows size is irrelevant. In this case the 
    `coeffs` parameter passed to the function would be:

        np.array([1, 0, 1, 1, 0.25])
        
    In this way the contribution of that parameter to the fitness will
    always be 0. Note that you may also want to set the corresponding 
    exponent to 0 as well. This may lead to a faster calculation. Notice
    that the carrot term has the coeffiecient set to 0.25. This gives
    it the same weighing as the other terms in the penalty term. The 
    difference is due to the fact that the penality term's exponents are
    summed first and then used in 1/x, while the carrot term is not 
    summed or passed through an inverse function.
  
    Parameters
    ----------
    macro_mol : Cage
        The cage whose fitness is to be calculated.
        
    population : Population
        The population in which `macro_mol` is held. Used for scaling 
        purposes.
        
    target_cavity : float
        The desried size of the cage's pore.
        
    macromodel_path : str
        The Schrodinger directory path.

    target_window : float (default = None)
        The desired radius of the largest window of the cage. If 
        ``None`` then `target_cavity` is used.
        
    coeffs : numpy.array (default = None)
        An array holding the coeffients A to N in equation (2).
        
    exponents : numpy.array (default = None)
        An array holding the exponents a to n in equation (2).
        
    energy_params : dict
        A dictionary holding the name arguments and values for the 
        ``Energy.pseudoformation()`` function.
    
    Returns
    -------
    The fitness function will be called on each cage twice. Once to get
    the unscaled values of all the var parameters. This will return an 
    array. The second time the average of the unscaled values is used 
    for scaling. The scaled values are then used to calculate the final 
    fitness as a float. This is returned after the second call. More
    details in ``calc_fitness``.
    
    float
        The fitness of `macro_mol`.
        
    numpy.array
        A numpy array holding unsacled values of var1 to varX. This is 
        returned during the scaling procedure described in 
        ``calc_fitness``.

    Raises
    ------
    ValueError
        When the average of the negavite bond formation energy is 
        positive.

    """

    # Calculate the mean values of the parameters across the population.
    # Only if they have not been calculated already however.
    if not population._fitness_cache:
        print(('\n\n\nCalculating mean values'
               ' of fitness parameters.\n\n\n'))
        # These lines calculate the unscaled parameters of all 
        # population members and return the means of those parameters.
        p_calc_vars = partial(_calc_fitness_vars(
                               target_cavity=target_cavity,
                               target_window=target_window,
                               energy_params=energy_params))
                               
        mean = population.mean(p_calc_vars)
        # The means are saved to the cache.
        population._fitness_cache['mean'] = mean

    # At this point all members of the population will have had their
    # unscaled fitness parameters calculated.

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

def _calc_fitness_vars(macro_mol, target_cavity, 
                       target_window, energy_params):
    """
    Calculates the unscaled fitness parameters for ``cage()``.

    Parameters
    ----------
    macro_mol : Cage
        The  Cage instance which needs to have its fitness calculated.
        
    target_cavity : float
        The desried size of the cage's cavity.
        
    target_window : float
        The desired radius of the largest window of the cage. If 
        ``None`` then `target_cavity` is used.
    
    Modifies
    --------
    macro_mol._unscaled_fitness_vars : numpy.array
        Creates this attribute and places the unscaled fitness
        parameters within.
    
    macro_mol._fitness_fail : bool
        If the calculation of one of the fitness parameters failed this
        is ``True``, else ``False``.
    
    Returns
    -------
    numpy.array
        An array holding the unscaled fitness parameters of `macro_mol`.
    
    """
     
    # If the parameters have already been calculated for this 
    # `macro_mol` do not recalculate them.
    if hasattr(macro_mol, '_unscaled_fitness_vars'):
        return macro_mol._unscaled_fitness_vars
                   
    cavity_diff = abs(target_cavity - macro_mol.topology.cavity_size())

    window_diff = (abs(target_window - 
                      max(macro_mol.topology.windows)) if 
                      macro_mol.topology.windows is not None else None)
                      
    asymmetry = (macro_mol.topology.window_difference() if
                 macro_mol.topology.window_difference() is not None 
                                                        else None) 

    e_per_bond = macro_mol.energy.pseudoformation(**energy_params)
    e_per_bond /= macro_mol.topology.bond_made

    if e_per_bond < 0:
        ne_per_bond = e_per_bond
        pe_per_bond = 0
    else:
        ne_per_bond = 0
        pe_per_bond = e_per_bond
    
    fitness_vars = [cavity_diff,
                    window_diff,
                    asymmetry,
                    ne_per_bond,
                    pe_per_bond]
                    
    macro_mol._fitness_fail = True if None in fitness_vars else False

    macro_mol._unscaled_fitness_vars = np.array([
                    cavity_diff,
                    (window_diff if window_diff is not None else 0),
                    (asymmetry if asymmetry is not None else 0),
                    ne_per_bond,
                    pe_per_bond])
                                               
    return macro_mol._unscaled_fitness_vars

def cage_target(macro_mol, target_mol_file, macromodel_path, 
                rotations=0):
    """
    Calculates the fitness of a cage / target complex.
    
    Parameters
    ----------
    macro_mol : Cage
        The cage which is to have its fitness calculated,

    target_mol_file : str
        The full path of the .mol file hodling the target molecule
        placed inside the cage.
        
    macromodel_path : str
        The Schrodinger directory path.

    rotations : int (default = 0)
        The number of times the target should be randomly rotated within 
        the cage cavity in order to find the most stable conformation.
        
    Returns
    -------
    float
        The fitness value of `macro_mol`.
    
    """
    
    # If the cage already has a fitness value, don't run the
    # calculation again.
    if macro_mol.fitness:
        print('Skipping {0}'.format(macro_mol.prist_mol_file))
        return macro_mol.fitness
    
    # The first time running the fitness function create an instance
    # of the target molecule as a ``StructUnit``. Due to caching,
    # running the initialization again on later attempts will not 
    # re-initialize.
    target = StructUnit(target_mol_file, minimal=True)

    # This function creates a new molecule holding both the target
    # and the cage centered at the origin. It then calculates the 
    # energy of this complex and compares it to the energies of the
    # molecules when separate. The more stable the complex relative
    # to the individuals the higher the fitness.
    
    # Create rdkit instances of the target in the cage for each
    # rotation.        
    rdkit_complexes = _generate_complexes(macro_mol, target, rotations+1)
    
    # Optimize the strcuture of the cage/target complexes.
    macromol_complexes = []        
    for i, complex_ in enumerate(rdkit_complexes):
        # In order to use the optimization functions, first the data is 
        # loaded into a ``MacroMolecule`` instance and its .mol file is 
        # written to the disk.
        mm_complex = MacroMolecule.__new__(MacroMolecule)
        mm_complex.prist_mol = complex_
        mm_complex.prist_mol_file = macro_mol.prist_mol_file.replace(
                            '.mol', '_COMPLEX_{0}.mol'.format(i))
        mm_complex.write_mol_file('prist')
        
        optimization.macromodel_opt(mm_complex, no_fix=True,
                       macromodel_path=macromodel_path)
        macromol_complexes.append(mm_complex)
    
    # Calculate the energy of the complex and compare to the
    # individual energies. If more than complex was made, use the
    # most stable version.
    energy_separate = macro_mol.energy + target.energy
    energy_diff =  min(macromol_complex.energy - energy_separate for 
                            macromol_complex in macromol_complexes)
    
                       
    raw_fitness = np.exp(energy_diff*1e-5) + 1
    if raw_fitness > 1e10:
        raw_fitness = 1e10
        
    return raw_fitness
   
def _generate_complexes(macro_mol, target, number=1):
    """
    Yields rdkit instances of cage / target complexes.
    
    If multiple complexes are returned, they will be different via a
    random rotation accross the x, y and z axes.
    
    Parameters
    ----------
    macro_mol : Cage
        The cage used to form the complex.
        
    target : StructUnit
        The target used to form the complex.
        
    number : int (default = 1)
        The number of complexes to be returned.
        
    Yields
    ------
    rdkit.Chem.rdchem.Mol
        An rdkit instance holding the cage / target complex. 
    
    """

    # First place both the target and cage at the origin.
    macro_mol.set_position('prist', [0,0,0])
    target.set_position('prist', [0,0,0])
    
    # Get the position matrix of the target molecule.        
    og_pos_mat = target.position_matrix('prist')
    
    # Carry out every rotation and yield a complex for each case.
    for i in range(number):
        rot_target = copy.deepcopy(target)
        
        rot1 = np.random.rand() * 2*np.pi
        rot2 = np.random.rand() * 2*np.pi
        rot3 = np.random.rand() * 2*np.pi
        
        rot_mat1 = rotation_matrix_arbitrary_axis(rot1, [1,0,0])
        rot_mat2 = rotation_matrix_arbitrary_axis(rot2, [0,1,0])
        rot_mat3 = rotation_matrix_arbitrary_axis(rot3, [0,0,1])
        
        new_pos_mat = np.dot(rot_mat1, og_pos_mat)
        new_pos_mat = np.dot(rot_mat2, new_pos_mat)
        new_pos_mat = np.dot(rot_mat3, new_pos_mat)
        
        rot_target.set_position_from_matrix('prist', new_pos_mat)
        
        yield chem.CombineMols(macro_mol.prist_mol, rot_target.prist_mol)
    


    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
