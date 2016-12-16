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

Fitness functions should return the `macro_mol` instance they took as
an argument. Within the fitness function itself the attribute 
`macro_mol.fitness` should have the fitness value calculated by the 
function placed into it. This ensures that the functions can be 
parallelized.

A scaling function can be defined in ``scaling.py`` if some sort of 
normalization of fitness values across the population is required. This 
is completely optional. Fitness functions should be atomic and depend
only on the individual whose fitness they are calculating.

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
import networkx as nx
import os
import multiprocessing as mp
import warnings
from collections import Counter

from .classes.exception import MacroMolError
from .classes.molecular import MacroMolecule, StructUnit, Energy
from .classes.topology.cage import CageTopology
from . import optimization
from .convenience_tools import (rotation_matrix_arbitrary_axis, 
                                matrix_centroid)

def _calc_fitness(func_data, population):
    """
    Calculates the fitness values of all members of a population.   
    
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
    p_func = partial(func, **func_data.params)

    # Apply the function to every member of the population, in parallel.
    with mp.get_context('spawn').Pool() as pool:
        evaluated = pool.map(p_func, population)

        # Make sure the cache is updated with the evaluated versions.
        for member in evaluated:
            member.update_cache()
        
    return evaluated

def _calc_fitness_serial(func_data, population):
    """
    Calculates the fitness values of all members of a population.  
    
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
            func(macro_mol, **func_data.params)
                
        except Exception as ex:
            MacroMolError(ex, macro_mol, 'During fitness calculation.')


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
    
    Modifies
    --------
    macro_mol.fitness : float
        Assigns a fitness to this attribute.
    
    Returns
    -------
    macro_mol
        The `macro_mol` with an integer between 0 (including) and 100 
        (excluding) as its fitness.

    """

    macro_mol.fitness = np.random.randint(1,10)
    return macro_mol    
    
# Provides labels for the progress plotter.
@_param_labels('Cavity Difference ','Window Difference ',
                'Asymmetry ', 'Positive Energy per Bond ', 
                'Negative Energy per Bond ')
def cage(macro_mol, target_cavity, target_window=None, 
         energy_params={'key':('rdkit', 'uff')}):
    """
    Calculates the fitness of a cage.
    
    This function is intended to be used with the normalization function
    ``carrots_and_sticks()`` defined in ``normalization.py``.

    The fitness function creates a tuple of 2 arrays. The first array
    holds parameters of `macro_mol` which contribute to a high fitness.
    The second array holds parameters of `macro_mol` which cause a low
    fitness.

    The parameter which indicates high fitness is the negative formation 
    energy per bond made. The remaining parameters indicate low fitness
    and are:
        1) `cavity_diff` - the difference between the cavity of 
           `macro_mol` and the `target_cavity`.
        2) `window_diff` - the difference between the largest window of
           `macro_mol` and `target_window`.
        3) `asymmetry` - the sum of the size differences of all the
           windows in `macro_mol`.
        4) `pos_eng_per_bond` - The postive formation energy of 
           `macro_mol` per bond made.
  
    Parameters
    ----------
    macro_mol : Cage
        The cage whose fitness is to be calculated.
        
    target_cavity : float
        The desried diameter of the cage's pore.

    target_window : float (default = None)
        The desired diameter of the largest window of the cage. If 
        ``None`` then `target_cavity` is used.
        
    energy_params : dict
        A dictionary holding the name arguments and values for the 
        ``Energy.pseudoformation()`` function.
    
    Modifies
    --------
    macro_mol.fitness_fail : bool
        The function sets this to ``True`` if one of the parameters 
        was not calculated. ``False`` if every parameter was calculated
        successfully.
        
    macro_mol.unscaled_fitness : tuple of 2 numpy.arrays
        The first numpy array holds the value of the negative energy
        per bond made. The second array holds the remaining parameters
        described above.   
        
    macro_mol.progress_params : list
        Places the calculated parameters in a single list. The order
        corresponds to the arguments in the ``_param_labels()`` 
        decorator applied to this function.
        
    Returns
    -------
    macro_mol : Cage
        The `macro_mol` with its fitness parameters calculated.

    """

    # Prevents warnings from getting printed when using multiprocessing.
    warnings.filterwarnings('ignore')
    try:
         
        # If the parameters have already been calculated for this 
        # `macro_mol` do not recalculate them.
        if macro_mol.unscaled_fitness:
            return macro_mol
                       
        if target_window is None:
            target_window = target_cavity                       
                       
        cavity_diff = abs(target_cavity - 
                          macro_mol.topology.cavity_size())
    
        if macro_mol.topology.windows is not None:
            window_diff = abs(target_window - 
                              max(macro_mol.topology.windows))
        else:
            window_diff  = None
            
        if  macro_mol.topology.window_difference() is not None:             
            asymmetry = macro_mol.topology.window_difference()
        else:
            asymmetry = None
    
        e_per_bond = macro_mol.energy.pseudoformation(**energy_params)
        e_per_bond /= macro_mol.topology.bonds_made
    
        if e_per_bond < 0:
            ne_per_bond = abs(e_per_bond)
            pe_per_bond = 0
        else:
            ne_per_bond = 0
            pe_per_bond = e_per_bond
        
        macro_mol.progress_params = [cavity_diff, window_diff, 
                                   asymmetry, pe_per_bond, -ne_per_bond]  
                        
        macro_mol.fitness_fail = (True if None in 
                                  macro_mol.progress_params else False)
    
        macro_mol.unscaled_fitness = (np.array([ne_per_bond]),
                        np.array([cavity_diff,
                        (window_diff if window_diff is not None else 0),
                        (asymmetry if asymmetry is not None else 0),
                        pe_per_bond]))
        
          
        
        return macro_mol

    except Exception as ex:
        MacroMolError(ex, macro_mol, "During fitness calculation")
        return macro_mol

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

@_param_labels('Negative Binding Energy', 'Positive Binding Energy', 
               'Asymmetry') 
def cage_c60(macro_mol, target_mol_file, 
             macromodel_path, n5fold, n2fold, min_cavity=None):
    """
    Calculates the fitness of a cage / C60 complex.
    
    The difference between this function and `cage_target()` is that
    the rotations are specifically aimed at sampling C60 entirely and
    systematically. Rather than the random sampling of the other
    function.
    
    Parameters
    ----------
    macro_mol : Cage
        The cage which is to have its fitness calculated,

    target_mol_file : str
        The full path of the .mol file hodling the target molecule
        placed inside the cage.
        
    macromodel_path : str
        The Schrodinger directory path.

    n5fold : int
        The number of rotations along the 5-fold axis of symmetry.
        
    n2fold : int
        The number of rotations along the 2 fold axis of symmetry per
        rotation along the 5-fold axis.

    Modifies
    --------
    macro_mol.fitness : float
        Places a fitness value into this attribute

    Returns
    -------
    macro_mol
        The `macro_mol` with its fitnes value calculated
    
    """
    
    warnings.filterwarnings('ignore')
    try:
        # If the cage already has a fitness value, don't run the
        # calculation again.
        if macro_mol.unscaled_fitness:
            print('Skipping {0}'.format(macro_mol.prist_mol_file))
            return macro_mol
     
        if min_cavity and min_cavity < macro_mol.topology.cavity_size():
            macro_mol.fitness = 1e-4
            return macro_mol
           
        # Make a copy version of `macro_mol` which is unoptimizted.
        unopt_macro_mol = copy.deepcopy(macro_mol)
        unopt_macro_mol.topology.final_sub()
        
        
        # The first time running the fitness function create an instance
        # of the target molecule as a ``StructUnit``. Due to caching,
        # running the initialization again on later attempts will not 
        # re-initialize.
        target = StructUnit(target_mol_file, minimal=True)
        _, molname = os.path.split(macro_mol.prist_mol_file)
        
        # Write a copy of the target for each macro_mol. So that parallel
        # energy calculations don't clash. There is room for optimization
        # here.
        molname, ext = os.path.splitext(molname)
        target.prist_mol_file = os.path.join(os.getcwd(), 
                                             molname+"target"+ext)    
        target.write_mol_file('prist')
        
    
        # This function creates a new molecule holding both the target
        # and the cage centered at the origin. It then calculates the 
        # energy of this complex and compares it to the energies of the
        # molecules when separate. The more stable the complex relative
        # to the individuals the higher the fitness.
        
        # Create rdkit instances of the target in the cage for each
        # rotation.        
        rdkit_complexes = _c60_rotations(unopt_macro_mol, target, 
                                         n5fold, n2fold)
    
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
            mm_complex.optimized = False
            mm_complex.energy = Energy(mm_complex)
            optimization.macromodel_opt(mm_complex, no_fix=True,
                           macromodel_path=macromodel_path)
            macromol_complexes.append(mm_complex)
    
        # Calculate the energy of the complex and compare to the
        # individual energies. If more than complex was made, use the
        # most stable version.
        energy_separate = (macro_mol.energy.macromodel(16, macromodel_path) + 
                            target.energy.macromodel(16, macromodel_path))
        min_eng_cmplx = min(macromol_complexes, 
                    key=lambda x : x.energy.macromodel(16, macromodel_path))                        
    
        binding_energy = (min_eng_cmplx.energy.values[('macromodel', 16)] - 
                           energy_separate)
            
        if binding_energy > 0:
            pos_be = binding_energy
            neg_be = 0
        else:
            pos_be = 0
            neg_be = abs(binding_energy)

        frag1, frag2 = chem.GetMolFrags(min_eng_cmplx.prist_mol, 
                                        asMols=True,
                                        sanitizeFrags=False)
                                      
        cage_counter = Counter(x.GetAtomicNum() for x in 
                                macro_mol.prist_mol.GetAtoms())
        frag_counters = [(frag1, Counter(x.GetAtomicNum() for x in 
                                frag1.GetAtoms())),

                        (frag2, Counter(x.GetAtomicNum() for x in 
                                frag2.GetAtoms()))]
    
        cmplx_cage_mol = next(frag for frag, counter in frag_counters if 
                            counter == cage_counter)
        
        cmplx_cage = MacroMolecule.__new__(MacroMolecule)
        cmplx_cage.prist_mol = cmplx_cage_mol
        cmplx_cage.topology = type(macro_mol.topology)(cmplx_cage)
        cmplx_cage.prist_mol_file = macro_mol.prist_mol_file.replace(
                         '.mol', '_COMPLEX_{0}_no_target.mol'.format(i))
        cmplx_cage.write_mol_file('prist')
        
    
        if cmplx_cage.topology.window_difference() is not None:             
            asymmetry = macro_mol.topology.window_difference()
        else:
            asymmetry = None        
    
        
        macro_mol.progress_params = [neg_be, pos_be, asymmetry]        
        
        macro_mol.fitness_fail = (True if None in 
                                   macro_mol.progress_params else False)
    
        macro_mol.unscaled_fitness = (
                              np.array([neg_be]),
                              np.array([pos_be, 
                              (asymmetry if asymmetry is not None else 0)]))
    
        return macro_mol
        
    except Exception as ex:
        MacroMolError(ex, macro_mol, "During fitness calculation.")
        return macro_mol
    
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
    
def _c60_rotations(macro_mol, c60, n5fold, n2fold):
    """
    Rotates C60 about its axes of symmetry and places it in `macro_mol`.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The cage which should have C60 placed inside it.
        
    c60 : StructUnit
        A StructUnit instance of C60.
        
    n5fold : int
        The number of rotations along the 5-fold axis of symmetry.
        
    n2fold : int
        The number of rotations along the 2 fold axis of symmetry per
        rotation along the 5-fold axis.
        
    Yields
    ------
    rdkit.Chem.rdchem.Mol
        An rdkit instance holding the cage / C60 complex. 
    
    """
    
    
    macro_mol.set_position('prist', [0,0,0])
    c60.set_position('prist', [0,0,0])
    
    # Step 1: Align the 5 membered ring with the z-axis.
    
    # Find a the ids of atoms in a membered ring.
    g = c60.graph('prist')
    ids = next(x for x in nx.cycle_basis(g) if len(x) == 5)
    # Place the coordinates of those atoms in a matrix.
    ring_matrix = np.matrix(
                        [c60.atom_coords('prist', id_) for id_ in ids])

    # Get the centroid of the ring.    
    ring_centroid = matrix_centroid(ring_matrix)
    # Align the centroid of the ring with the z-axis.
    c60.set_orientation('prist', ring_centroid, [0,0,1])
    aligned_c60 = copy.deepcopy(c60)
    
    # Step 2: Get the rotation angles and apply the rotations. Yield 
    # the resulting complex.
    
    # Get the angles of the 5 and 2 fold rotations.
    angles5fold = np.arange(0, 72/180*np.pi, 72/180*np.pi/n5fold)
    angles2fold = np.arange(0, np.pi, np.pi/n2fold)
    
    for angle5 in angles5fold:
        for angle2 in angles2fold:
            buckyball = copy.deepcopy(aligned_c60)
            buckyball.rotate('prist', angle5, [0,0,1])
            buckyball.rotate('prist', angle2, [0,1,0])
            yield chem.CombineMols(macro_mol.prist_mol, buckyball.prist_mol)

    
    
    

    
    
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
