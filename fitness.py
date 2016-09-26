import numpy as np
import rdkit.Chem as chem
import itertools as it
import copy

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
            macro_mol.fitness = 0
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

def cage(macro_mol, target_size, coeffs=None, exponents=None):
    if macro_mol.fitness:
        print('Skipping {0}'.format(macro_mol.prist_mol_file))
        return macro_mol.fitness
    
    if macro_mol.topology.windows is None:
        return 0

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

    return 1/np.sum(fitness_value)   
    
def cage_target(cage, target_mol_file, *, 
                macromodel_path, rotate=False):
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
    
                       
    return np.exp(energy_diff) + 1
   
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
    


    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    