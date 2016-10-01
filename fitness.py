import numpy as np
import rdkit.Chem as chem
import itertools as it
import copy
from scipy.stats import logistic
import os
from inspect import signature

from .classes.exception import MacroMolError
from .classes.molecular import MacroMolecule, StructUnit
from .optimization import macromodel_opt
from .convenience_functions import rotation_matrix_arbitrary_axis

def calc_fitness(func_data, population):
    """
    Calculates the fitness values of all members of a population.
    
    A fitness function should take a ``MacroMolecule`` instance and
    return a number representing its fitness. The assignement to the
    `fitness` attribute of a population member happens here, not by the
    fitness function.    
    
    Extending MMEA: Adding fitness functions
    ----------------------------------------
    To add a new fitness function simply write it as a function in this
    module. It will need to take the ``MacroMolecule`` instance as its
    first argument. 
    
    Some fitness functions will need to have certain internal values 
    normalized. In this section, the ``cage`` fitness function is used
    as an example.
    
    The cage fitness function has the form
    
        (1) fitness = A*(var1^a) + B(var2^b) + C(var3^c) + D(var4^d).
        
    Here var1 to var4 represent the parameters of a cage which factor
    into its fitness. If raw values of these parameters are used, due
    to different units and naturally different orders of magnitude one
    parameter may have an undue influence on the fitness only because
    of its units. The goal is for fitness function contributions to be
    determined soley by the coeffiecients A to D and exponents a to d.
    
    In order to do this the average value of var1 throughout the entire
    population is calculated first and then used as a scaling factor.
    The variables in the fitness function therefore have the form
    
        (2) var = var_ind / <var>,
        
    where var can be any of var1 to var4 in equation (1), var_ind is the
    variable of for that individual and <var> is the average of that 
    variable throughout the population.
    
    In order to calculate the fitness in this way the fitness function
    needs to be applied twice to each member of the population. The
    first loop calculates var_ind and <var> and the second loop 
    calculates var. See the implementation.
    
    Fitness functions which require this scaling procedure will need to
    have a keyword argument called `means` which is default initialized
    to ``None``. In order to make use of this they will also have to
    have the general form of equation (1). Nothing else is reqiured
    but they should allow the user to supply array holding the c
    oefficients and exponents. See the implementation of ``cage``.
    
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
    
    use_means = 'means' in signature(func).parameters.keys() 
    
    if use_means:

        var_sum = 0
        for ind in population:
            try:
                unscaled = func(ind, **func_data.params)
            except:
                unscaled = 0
            var_sum = np.add(unscaled, var_sum)

        var_avg = np.divide(var_sum, len(population))

    print(var_avg)

    # Apply the function to every member of the population.
    for macro_mol in population:
        try:
            if use_means:
                macro_mol.fitness = func(macro_mol, means=var_avg,
                                         **func_data.params)
            else:
                macro_mol.fitness = func(macro_mol, **func_data.params)
                
        except Exception as ex:
            MacroMolError(ex, macro_mol, 'During fitness calculation.')
            macro_mol.topology.windows = None

        print(macro_mol.fitness, '-', macro_mol.prist_mol_file)            
            
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

def cage(macro_mol, target_size, 
         coeffs=None, exponents=None, means=None):
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
        
        1) `cavity_diff` - the difference between the cage's cavity and
           the desired cavity size
        2) `window_area_diff` - the difference between the size of the 
           window of the cage and the size required to allow the target
           to enter the cavity
        3) `asymmetry` - sum of the difference between the size of of
           windows of the same type
        4) `neg_eng_per_bond` - the energy of the cage divided by the 
           number of bonds for during assembly, when the total energy
           of the cage is < 0. This is a measure of the stability or 
           strain of a cage.
        5) Same as 4) but used when total energy of the cage is > 0.
        
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
    the `window_area_diff` parameter. This may be because cages which 
    form via templating are also to be considered. In this case the 
    `coeffs` parameter passed to the function would be:

        np.array([1,0,1,1, 0.25])
        
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
        
    target_size : float
        The desried size of the cage's pore.
        
    coeffs : numpy.array (default = None)
        An array holding the coeffients A to N in equation (2).
        
    exponents : numpy.array (default = None)
        An array holding the exponents a to n in equation (2).
        
    means : numpy.array (default = None)
        A numpy array holding the mean values of var1 to varx over the
        populatoin. Used in the scaling procedure decribed in 
        ``calc_fitnes``.
    
    Returns
    -------
    float
        The fitness of `macro_mol`.
        
    numpy.array
        A numpy array holding values of var1 to varx. This is returned
        during the scaling procedure described in ``calc_fitness``.

    """
    
    if macro_mol.fitness:
        print('Skipping {0}'.format(macro_mol.prist_mol_file))
        return macro_mol.fitness
    
    if macro_mol.topology.windows is None:
        return 1

    if means is not None:
        if coeffs is None:
            coeffs = np.array([1,1,1,1,0.2])
            
        if exponents is None:
            exponents = np.array([1,1,1,1,1])  
        
        # Make sure you are not dividing by 0.
        for i, x in enumerate(means):
            if x == 0:
                means[i] = 1
        
        scaled = np.divide(macro_mol.unscaled_fitness_vars, means)
        # Delete the attribute once is no longer necessary.
        delattr(macro_mol, 'unscaled_fitness_vars')      
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

    cavity_diff = abs(target_size - macro_mol.topology.cavity_size())

    target_window_area = np.square(target_size)
    window_area = np.square(max(macro_mol.topology.windows))
    window_area_diff = abs(target_window_area - window_area)    
      
    asymmetry = macro_mol.topology.window_difference(500)
    
    energy_per_bond = macro_mol.energy / macro_mol.topology.bonds_made
    if energy_per_bond < 0:
        neg_eng_per_bond = energy_per_bond
        pos_eng_per_bond = 0
    else:
        neg_eng_per_bond = 0
        pos_eng_per_bond = energy_per_bond

    unscaled =  np.array([
                     cavity_diff, 
                     window_area_diff,                                                          
                     asymmetry,
                     pos_eng_per_bond,
                     neg_eng_per_bond
                     ])
    macro_mol.unscaled_fitness_vars = unscaled
    return unscaled
    
def cage_target(cage, target_mol_file, target_size, *, macromodel_path, 
                rotate=False, min_window_size=0):
    """
    Calculates the fitness of a cage / target complex.
    
    Depending of `rotate` a different number of cage / target 
    complexes will be generated. When `rotate` is ``True``, the 
    target molecule is rotated along the x, y and z axes and each
    rotation forms a new complex. The lowest energy complex is used
    for the fitness calculation. When `rotate` is ``False`` not 
    rotation takes place.
    
    To see which rotations take place see the documentation of the
    `generate_complexes` method.
    
    Parameters
    ----------
    cage : Cage
        The cage which is to have its fitness calculated,

    target_mol_file : str
        The full path of the .mol file hodling the target molecule
        placed inside the cage.
        
    target_size : float
        The minimum size which the cage cavity needs to be in order to
        encapsulate the target.
        
    min_window_size : float (default = 0)
        The smallest windows size allowing the target to enter the cage.
        Default is 0, which implies that there is no minimum. This can
        occur when the target acts a template for cage assembly.

    rotate : bool (default = False)
        When ``True`` the target molecule will be rotated inside the
        cavity of the cage. When ``False`` only the orientation
        within the .mol file is used.
        
    macromodel_path : str (keyword-only)
        The Schrodinger directory path.
        
    Returns
    -------
    float
        The fitness value of `cage`.
    
    """
    
    # If the cage already has a fitness value, don't run the
    # calculation again.
    if cage.fitness:
        print('Skipping {0}'.format(cage.prist_mol_file))
        return cage.fitness
            
    # If the size of the window is too small for the target molecule to
    # enter the cage return the minimum fitness value.
    if max(cage.topology.windows) < min_window_size:
        return 1
        
    # If the size of the cage cavity is too small for the target, return
    # the minimum fitness value.
    if cage.topology.cavity_size() < target_size:
        return 1
    
    # The first time running the fitness function create an instance
    # of the target molecule as a ``StructUnit``. Due to caching,
    # running the initialization again on later attempts will not 
    # re-initialize. Normally ``StructUnit`` instances do not have
    # the attribute `optimized`, as these structures do not need to
    # be optimized. However, in this case it is necessary to know
    # the energy of the target molecule as part of the fitness 
    # calculation. Giving ``StructUnit`` the `optimized` attribute
    # fools the optimization function into running on ``StructUnit``
    # despite it being designed for ``MacroMolecule`` instances.
    target = StructUnit(target_mol_file, minimal=True)
    if not hasattr(target, 'optimized'):
        target.optimized = False
        macromodel_opt(target, no_fix=True, 
                       macromodel_path=macromodel_path)

    # This function creates a new molecule holding both the target
    # and the cage centered at the origin. It then calculates the 
    # energy of this complex and compares it to the energies of the
    # molecules when separate. The more stable the complex relative
    # to the individuals the higher the fitness.
    
    # Create rdkit instances of the target in the cage for each
    # rotation.        
    rdkit_complexes = list(_generate_complexes(cage, target, rotate))
    
    # Place the rdkit complexes into a new .mol file and use that 
    # to make a new ``StructUnit`` instance of the complex.
    # ``StructUnit`` is the class of choice here because it can be
    # initialized from a .mol file. ``MacroMolecule`` requires
    # things like topology and building blocks for initialization.
    macromol_complexes = []        
    for i, complex_ in enumerate(rdkit_complexes):
        # First the data is loaded into a ``MacroMolecule`` instance
        # as this is the class which is able to write rdkit
        # instances to .mol files. Note that this ``MacroMolecule``
        # instance is just a dummy and only holds the bare minimum
        # information required to write the complex to a .mol file.
        mm_complex = MacroMolecule.__new__(MacroMolecule)
        mm_complex.prist_mol = complex_
        mm_complex.prist_mol_file = cage.prist_mol_file.replace(
                            '.mol', '_COMPLEX_{0}.mol'.format(i))
        mm_complex.write_mol_file('prist')
        
        # Once the .mol file is written load it into a 
        # ``StructUnit`` instance.
        macromol_complex = StructUnit(mm_complex.prist_mol_file, 
                                      minimal=True)
        # In order to optimize a ``StructUnit`` instance the
        # `optimized` attribute needs to be added.
        macromol_complex.optimized = False
        macromodel_opt(macromol_complex, no_fix=True,
                       macromodel_path=macromodel_path)
        macromol_complexes.append(macromol_complex)
    
    # Calculate the energy of the complex and compare to the
    # individual energies. If more than complex was made, use the
    # most stable version.
    energy_separate = cage.energy + target.energy
    energy_diff =  min(energy_separate - macromol_complex.energy for 
                            macromol_complex in macromol_complexes)
    
                       
    raw_fitness = np.exp(energy_diff*1e-5) + 1
    if raw_fitness > 1e10:
        raw_fitness = 1e10
        
    return raw_fitness
   
def _generate_complexes(cage, target, rotate):
    """
    Yields rdkit instances of cage / target complexes.
    
    Parameters
    ----------
    cage : Cage
        The cage used to form the complex.
        
    target : StructUnit
        The target used to form the complex.
        
    rotate : bool
        When ``True`` the target molecule will undergo rotations
        within the cage cavity and a complex will be yielded for 
        each configuration.
        
    Yields
    ------
    rdkit.Chem.rdchem.Mol
        An rdkit instance holding the cage / target complex. 
    
    """
    
    # Define the rotations which are to be used on the target.
    if rotate:        
        rotations = [0, np.pi/2, np.pi, 3*np.pi/2]
    else:
        rotations = [0]

    # First place both the target and cage at the origin.
    cage.set_position('prist', [0,0,0])
    target.set_position('prist', [0,0,0])
    
    # Get the position matrix of the target molecule.        
    og_pos_mat = target.position_matrix('prist')
    
    # Carry out every rotation and yield a complex for each case.
    for rot1, rot2, rot3 in it.combinations_with_replacement(
                                                    rotations, 3):
        rot_target = copy.deepcopy(target)
        rot_mat1 = rotation_matrix_arbitrary_axis(rot1, [1,0,0])
        rot_mat2 = rotation_matrix_arbitrary_axis(rot2, [0,1,0])
        rot_mat3 = rotation_matrix_arbitrary_axis(rot3, [0,0,1])
        
        new_pos_mat = np.dot(rot_mat1, og_pos_mat)
        new_pos_mat = np.dot(rot_mat2, new_pos_mat)
        new_pos_mat = np.dot(rot_mat3, new_pos_mat)
        
        rot_target.set_position_from_matrix('prist', new_pos_mat)
        
        yield chem.CombineMols(cage.prist_mol, rot_target.prist_mol)
    


    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    