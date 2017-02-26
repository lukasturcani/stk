"""
Module for defining fitness functions.

A note on how fitness values are calculated.
--------------------------------------------
Calculation of fitness values can be a multi-step process. The goal is
to place a numerical value between 0 (exclusive) and infinity in the
`fitness` attribute of a MacroMolecule instance. A thing to keep in
mind is that once a value is placed into the `unscaled_fitness`
attribute it never changes (for a given fitness function and molecule).
However, values placed into the `fitness` attribute are subject to
change. This is because the outcome of a normalization procedure can
depend on other members of the population, not just the one being
evaulated. As a result this can be different each generation.

There are 2 ways in which calculation of a `fitness` value can be
achieved. First, using only a fitness function. Second, using a fitness
function and normalization functions.

The case when only a fitness function is used is simple. All fitness
functions take a MacroMoleule instance as an argument and return the
value of its fitness. MMEA then automatically puts this returned value
into the `unscaled_fitness` attribute. Next MMEA copies this value into
the `fitness` attribute.

So what happens if normalization functions are used?

The first thing to note is that multiple normalization functions can
be applied sequentially. Each normalization function replaces the
previous value in the `fitness` attribute. Normalization functions do
not manipulate or intarct with the `unscaled_fitness` attribute in any
way. Before the first normalization function is applied, MMEA
automatically copies the value in `unscaled_fitness` into `fitness`.

The normalization functions then place various values into the
`fitness` attribute. Only the last normalization function needs to
place a value between 0 (exclusive) and infinity in the `fitness`
attribute. This is useful if the fitness function calculates the values
a number of fitness parameters (such as energy, molecular weight, etc.)
and then the normalization functions combine them into a single number.

Each generation, before all the normalization functions are reapplied,
MMEA automatically copies the value in `unscaled_fitness` into
`fitness`. Then the sequence of normalization functions is applied
again.

While fitness functions are only applied once per molecule,
normalization functions are reapplied each generation.

Extending MMEA: Adding fitness functions.
-----------------------------------------

Example: random_fitness()

To add a new fitness function simply write it as a function in this
module. It will need to take the ``MacroMolecule`` instance as its
first argument and this argument should be called ``macro_mol``. The
purpose of this is to help users identify which arguments are handled
automatically by MMEA and which they need to define in the input file.
The convention is that if the fitness function takes an argument called
``macro_mol`` they do not have to specify that argument in the input
file.

A fitness function must return the value which holds the fitness of the
molecule taken as an argument. If a fitness function is meant to be
paired with a normalization funtion it can return any value or object
it likes. Just as long as the normalization functions know how to deal
with it and convert it to a number.

A fitness function may be complex and may not fit neatly into a single
function. For example, the ``cage_target()`` fitness function needs to
call ``_generate_complexes()`` in order to sample various conformations
before outputting a fitness value. This is fine. Define helper
functions such as ``_generate_complexes()`` within this module but make
sure they are private. This means that names of helper functions begin
with a leading underscore.

A note on plotting.
-------------------
As mentioned before some fitness functions may be complex and as a
result manipulate all sorts of data. Typically, in order to measure the
progress of a GA, the fitness values in the population are tracked
across generations. However, let's say that some hypothetical fitness
function also calculates the energies of molecules. It may be quite
interesting plot the evolution of energies across generations too. If
this is the case the fitness function may assign to the
`progress_params` attribute of `macro_mol`:

    macro_mol.progress_params = [mol_energy]

Now a plot showing the change in `mol_energy` across generations will
be made too, along with the plot showing the changes in fitness.

What if two things are needed to be kept track of?

    macro_mol.progress_params = [mol_energy, mol_radius]

Great, now a progress plot for each of the variables will be made.

How will the y axes be labelled in each plot?
The decorator `_param_labels()` exists for this.

Let's create a basic outline of a some fitness function:

    @_param_labels('Molecule Energy / J mol-1', 'Mol Radius / m-9')
    def this_is_the_fitness_function(macro_mol, some_param):
        ...
        calculate_stuff()
        ...
        macro_mol.progress_params = [mol_energy, mol_radius]
        ...
        return fitness_value

If this function is used in the GA, a progress plot will be made for
each of the `progress_params` and they will have their y-axes labelled
'Molecule Energy / J mol-1' and 'Molecule Radius / m-9', respectively.

"""

import numpy as np
import rdkit.Chem as chem
import copy
from functools import partial, wraps
import networkx as nx
import multiprocessing as mp
import warnings
from collections import Counter
import os
from os.path import join
from uuid import uuid4

from ..convenience_tools import (matrix_centroid,
                                 FunctionData, MolError,
                                 rotation_matrix_arbitrary_axis)

from ..molecular import (MacroMolecule,
                         StructUnit, Energy, optimization,
                         func_key)

def _calc_fitness(func_data, population):
    """
    Calculates the fitness values of all members of a population.

    Parameters
    ----------
    func_data : FunctionData
        A ``FunctionData`` instance representing the chosen fitness
        function and any additional parameters it may require.

    population : Population
        The population whose members must have their fitness
        calculated.

    Returns
    -------
    list
        The members of `population` which have had their fitness
        calculated.

    """

    # Get the fitness function object.
    func = globals()[func_data.name]
    # Make sure it won't raise errors while using multiprocessing.
    p_func = _FitnessFunc(partial(func, **func_data.params))

    # Apply the function to every member of the population, in
    # parallel.
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
        The population whose members must have their fitness
        calculated.

    Returns
    -------
    None : NoneType

    """

    # Get the fitness function object.
    func = globals()[func_data.name]


    # Apply the function to every member of the population.
    for macro_mol in population:
        _FitnessFunc(func(macro_mol, **func_data.params))


def _param_labels(*labels):
    """
    Adds the `param_labels` attribute to a fitness function.

    The point of this decorator is described in the module level
    docstring.

    Parameters
    ----------
    labels : tuple
        List of strings about the fitness labels used for plotting
        EPPs. The order of the strings should represent the order of
        the fitness ``vars`` in the fitness funciton. In practice it
        should correspond to the order of the ``coeffs`` or
        ``exponents`` parameters given to the fitness function.

    Returns
    -------
    func
        Decorated function.

    """

    def add_labels(func):
        func.param_labels = labels
        return func

    return add_labels

class _FitnessFunc:
    """
    A decorator for fitness functions.

    This decorator is applied to all fitness functions automatically in
    _calc_fitness(). It should not be applied explicitly when defining
    the functions.

    The decorator prevents fitness functions from raising if
    they fail (necessary for multiprocessing), prevents them from
    being run twice on the same molecule and stores the value returned
    by them in the `unscaled_fitness` dictionary.

    """

    def __init__(self, func):
        wraps(func)(self)

    def __call__(self, macro_mol, *args,  **kwargs):
        func_name = self.__wrapped__.func.__name__
        try:
            if (macro_mol.failed or
                func_name in macro_mol.unscaled_fitness):

                print('Skipping {0}'.format(macro_mol.name))
                return macro_mol

            val = self.__wrapped__(macro_mol, *args, **kwargs)
            macro_mol.unscaled_fitness[func_name] = val
            return macro_mol

        except Exception as ex:
            macro_mol.failed = True
            macro_mol.unscaled_fitness[func_name] = None
            MolError(ex, macro_mol, "During fitness calculation")
            return macro_mol

def random_fitness(macro_mol):
    """
    Returns a random fitness value.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to which a fitness value is to be assigned.

    Returns
    -------
    float
        A random postive number.

    """

    return abs(np.random.normal(50,20))

@_param_labels('var1', 'var2', 'var3', 'var4')
def random_fitness_vector(macro_mol):
    """
    Returns an array of random numbers.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to have its fitness calculated.

    Modifies
    --------
    macro_mol.progress_params : list
        The random numbers are also placed into this attribute.

    Returns
    -------
    numpy.array
        An array holding random numbers.

    """

    # Make a random fitness vector.
    f = abs(np.random.normal(50,20,4))
    # This multiplication ensures that the elements of the fitness
    # vector all have different oraders of magnitude and that some
    # are negative.
    f = np.multiply(f, np.array([0.01, 1, 10, -100]))
    macro_mol.progress_params = f.tolist()
    return f

def raiser(macro_mol, param1, param2=2):
    """
    Doens't calculate a fitness value, raises an error instead.

    This function is used to test that when fitness functions raise
    errors during multiprocessing, they are handled correctly.

    Parameters
    ---------
    param1 : object
        Dummy parameter, does nothing.

    param2 : object (default = 2)
        Dummy keyword parameter, does nothing.

    Returns
    -------
    This function does not return. It only raises.

    Raises
    ------
    Exception
        An exception is always raised.

    """

    raise Exception('Raiser fitness function used.')


# Provides labels for the progress plotter.
@_param_labels('Cavity Difference ','Window Difference ',
                'Asymmetry ', 'Energy per Bond ')
def cage(macro_mol, target_cavity, target_window=None,
         pseudoformation_params=
         { 'func' : FunctionData('rdkit', forcefield='mmff') }):
    """
    Returns the fitness vector of a cage.

    The fitness vector consists of the following properties in the
    listed order

        1) `cavity_diff` - the difference between the cavity of
           `macro_mol` and the `target_cavity`.
        2) `window_diff` - the difference between the largest window of
           `macro_mol` and `target_window`.
        3) `asymmetry` - the sum of the size differences of all the
           windows in `macro_mol`.
        4) `eng_per_bond` - The formation energy of `macro_mol` per
           bond made.

    Parameters
    ----------
    macro_mol : Cage
        The cage whose fitness is to be calculated.

    target_cavity : float
        The desired diameter of the cage's pore.

    target_window : float (default = None)
        The desired diameter of the largest window of the cage. If
        ``None`` then `target_cavity` is used.

    pseudoformation_params : dict (default =
          { 'energy_func' : FunctionData('rdkit', forcefield='uff') })

        This fitness function calculates the formation energy using the
        ``Energy.pseudoformation()`` method. This parameter defines the
        arguments passed to this method via a dictionary. The name of
        the argument is the key and the value of the argument is the
        value.

        Default initialized arguments of Energy.pseudoformation() only
        need to be specified in `energy_params` if the user wishes to
        change the default value.

        To see what arguments the `Energy.pseudoformation()` method
        requires, try using the  `-h` option:

            python -m mmea -h energy

    Modifies
    --------
    macro_mol.failed : bool
        The function sets this to ``True`` if one of the parameters
        was not calculated.

    macro_mol.progress_params : list
        Places the calculated parameters in the list. The order
        corresponds to the arguments in the ``_param_labels()``
        decorator applied to this function.

    Returns
    -------
    numpy.array
        The numpy array holds the fitness vector described in this
        docstring.

    None : NoneType
        Returned if any fitness parameter failed to calculate.

    """

    # Prevents warnings from getting printed when using
    # multiprocessing.
    warnings.filterwarnings('ignore')

    if target_window is None:
        target_window = target_cavity

    cavity_diff = abs(target_cavity -
                      macro_mol.cavity_size())

    if macro_mol.windows is not None:
        window_diff = abs(target_window -
                          max(macro_mol.windows))
    else:
        window_diff  = None

    asymmetry = macro_mol.window_difference()

    print('\n\nCalculating complex energies.\n')
    e_per_bond = macro_mol.energy.pseudoformation(
                                           **pseudoformation_params)
    e_per_bond /= macro_mol.bonds_made

    macro_mol.progress_params = [cavity_diff, window_diff,
                                 asymmetry, e_per_bond]

    if None in macro_mol.progress_params:
        macro_mol.failed = True
        return None

    return np.array([cavity_diff, window_diff, asymmetry, e_per_bond])

@_param_labels('Binding Energy', 'Asymmetry')
def cage_target(macro_mol, target_mol_file,
                eng_func, opt_func, rotations=0):
    """
    Returns the fitness vector of a cage / target complex.

    The target is randomly rotated inside the cage's cavity and the
    most stable conformation found is used.

    The function calculates a fitness vector. The fitness vector
    consists of the following properties in the listed order,

            1) binding energy
            2) asymmetry

    Parameters
    ----------
    macro_mol : Cage
        The cage which is to have its fitness calculated,

    target_mol_file : str
        The full path of the .mol file hodling the target molecule
        placed inside the cage.

    energy_func : FunctionData
        A FunctionData object representing the energy function used to
        calculate energies.

    opt_func : FunctionData
        A FunctionData object representing the optimization function to
        be run on the generated complexes.

    rotations : int (default = 0)
        The number of times the target should be randomly rotated
        within the cage cavity in order to find the most stable
        conformation.

    Modifies
    --------
    macro_mol.progress_params : list
        Places the calculated parameters in the list. The order
        corresponds to the arguments in the ``_param_labels()``
        decorator applied to this function.

    macro_mol.failed : bool
        The function sets this to ``True`` if one of the parameters
        was not calculated.

    Returns
    -------
    numpy.array
        The numpy array holds the fitness vector described in this
        docstring.

    None : NoneType
        Returned if any fitness parameter failed to calculate.

    """

    return _cage_target(macro_mol,
                        target_mol_file, energy_func, opt_func,
                        FunctionData('_generate_complexes',
                                     rotations+1))

@_param_labels('Binding Energy', 'Asymmetry')
def cage_c60(macro_mol, target_mol_file,
             energy_func, opt_func, n5fold, n2fold):
    """
    Calculates the fitness vector of a cage / C60 complex.

    The difference between this function and `cage_target()` is that
    the rotations are specifically aimed at sampling C60 entirely and
    systematically. Rather than the random sampling of the other
    function.

    The fitness vector consists of the following properties in the
    listed order,

        1) binding energy
        2) asymmetry

    Parameters
    ----------
    macro_mol : Cage
        The cage which is to have its fitness calculated.

    target_mol_file : str
        The full path of the .mol file hodling the target molecule
        placed inside the cage.

    energy_func : FunctionData
        A FunctionData object representing the energy function used to
        calculate energies.

    opt_func : FunctionData
        A FunctionData object representing the optimization function to
        be run on the generated complexes.

    n5fold : int
        The number of rotations along the 5-fold axis of symmetry.

    n2fold : int
        The number of rotations along the 2 fold axis of symmetry per
        rotation along the 5-fold axis.

    Modifies
    --------
    macro_mol.progress_params : list
        Places the calculated parameters in the list. The order
        corresponds to the arguments in the ``_param_labels()``
        decorator applied to this function.

    macro_mol.failed : bool
        The function sets this to ``True`` if one of the parameters
        was not calculated.

    Returns
    -------
    numpy.array
        The numpy array holds the fitness vector described in this
        docstring.

    None : NoneType
        Returned if any fitness parameter failed to calculate.

    """
    return _cage_target(macro_mol,
                        target_mol_file, energy_func, opt_func,
                        FunctionData('_c60_rotations',
                                     n5fold=n5fold,
                                     n2fold=n2fold))

def _cage_target(macro_mol, target_mol_file,
                 energy_func, opt_func, rotation_func):
    """
    A general fitness function for calculating fitness of complexes.

    This function should be inherited by other fitness functions which
    define their own rotation function. For example ``cage_c60()`` and
    ``cage_target()``.

    The function returns a fitness vector consisting of the
    binding energy and asymmetry.

    Parameters
    ----------
    macro_mol : Cage
        The cage which is to have its fitness calculated.

    target_mol_file : str
        The full path of the .mol file hodling the target molecule
        placed inside the cage.

    energy_func : FunctionData
        A FunctionData object representing the energy function used to
        calculate energies.

    opt_func : FunctionData
        A FunctionData object representing the optimization function to
        be run on the generated complexes.

    rotation_func : FunctionData
        A FunctionData object representing the rotation function to be
        used.

    Modifies
    --------
    macro_mol.progress_params : list
        Places the calculated parameters in the list. The order
        corresponds to the arguments in the ``_param_labels()``
        decorator applied to this function.

    macro_mol.failed : bool
        The function sets this to ``True`` if one of the parameters
        was not calculated.

    Returns
    -------
    numpy.array
        The numpy array holds the fitness vector described in this
        docstring.

    None : NoneType
        Returned if any fitness parameter failed to calculate.

    """

    warnings.filterwarnings('ignore')

    # Create a folder to hold all the complexes generated by this
    # function.
    folder_path = _make_cage_target_folder()
    # Transform the FunctionData instances into functions.
    efunc = getattr(Energy, energy_func.name)
    ofunc = getattr(optimization, opt_func.name)
    rot_func = globals()[rotation_func.name]

    # Make a template name for all the complexes generated by this
    # function. If `macro_mol` has a `name` attribute use that as the
    # template. If `macro_mol` does not have a `name` the template will
    # just be a random, unique number.
    name = getattr(macro_mol, 'name', str(uuid4().int))

    # Make a copy version of `macro_mol` which is unoptimizted.
    unopt_macro_mol = copy.deepcopy(macro_mol)
    unopt_macro_mol.topology.build(unopt_macro_mol)

    # Create an instance of the target molecule as a ``StructUnit``.
    target = StructUnit(target_mol_file)

    # This function creates a new molecule holding both the target
    # and the cage centered at the origin. It then calculates the
    # energy of this complex and compares it to the energies of the
    # molecules when separate. The more stable the complex relative
    # to the individuals the higher the fitness.

    # Create rdkit instances of the target in the cage for each
    # rotation.
    rdkit_complexes = rot_func(macro_mol, target,
                               **rotation_func.params)

    # Optimize the strcuture of the cage/target complexes.
    macromol_complexes = []
    print('\n\nOptimizing complex structures.\n')
    for i, complex_ in enumerate(rdkit_complexes):
        # In order to use the optimization functions, first the data
        # is loaded into a ``MacroMolecule`` instance and its .mol
        # file is written to the disk.
        mm_complex = MacroMolecule.__new__(MacroMolecule)
        mm_complex.mol = complex_
        mm_complex.name = name + '_COMPLEX_{0}'.format(i)
        mm_complex.optimized = False
        mm_complex.energy = Energy(mm_complex)
        mm_complex.topology = macro_mol.topology
        mm_complex.building_blocks = macro_mol.building_blocks
        ofunc(mm_complex, **opt_func.params)
        macromol_complexes.append(mm_complex)

    # Calculate the energy of the complex and compare to the
    # individual energies. If more than complex was made, use the
    # most stable version.
    energy_separate = (efunc(macro_mol, **energy_func.params) +
                       efunc(target, **energy_func.params))

    print('\n\nCalculating complex energies.\n')
    min_eng_cmplx = min(macromol_complexes,
                        key=lambda x : efunc(x, **energy_func.params))

    # Write the most stable complex to a file.
    min_eng_cmplx.write(join(folder_path, min_eng_cmplx.name+'.mol'))

    ekey = func_key(efunc, min_eng_complx, energy_func.params)
    binding_energy = (min_eng_cmplx.energy.values[ekey] -
                                                    energy_separate)

    frag1, frag2 = chem.GetMolFrags(min_eng_cmplx.mol,
                                    asMols=True,
                                    sanitizeFrags=False)

    cage_counter = Counter(x.GetAtomicNum() for x in
                            macro_mol.mol.GetAtoms())
    frag_counters = [(frag1, Counter(x.GetAtomicNum() for x in
                            frag1.GetAtoms())),

                    (frag2, Counter(x.GetAtomicNum() for x in
                            frag2.GetAtoms()))]

    cmplx_cage_mol = next(frag for frag, counter in frag_counters if
                        counter == cage_counter)

    cmplx_cage = MacroMolecule.__new__(MacroMolecule)
    cmplx_cage.mol = cmplx_cage_mol
    cmplx_cage.topology = macro_mol.topology
    cmplx_cage.name = min_eng_cmplx.name + '_no_target'

    # Write the cage without the target to a file.
    cmplx_cage.write(join(folder_path, cmplx_cage.name+'.mol' ))

    asymmetry = macro_mol.window_difference()

    macro_mol.progress_params = [binding_energy, asymmetry]


    if None in macro_mol.progress_params:
        macro_mol.failed = True
        return None

    return np.array([binding_energy, asymmetry])

def _make_cage_target_folder():
    """
    Creates a folder to store molecules made by _cage_target().

    The function creates a folder called `cage_target`.
    Inside will be any complexes formed by the function _cage_target().
    The folder will be placed in the current working directory, or
    if the GA is running, 1 above the current working dirctory. This
    prevents the generated molecules from being cleaned up by the GA
    when it's finished.

    Returns
    -------
    str
        The path of the ``cage_target`` folder.

    """

    dir_path = os.getcwd()
    if join('output', 'scratch') in dir_path:
        dir_path = dir_path.replace('scratch', 'cage_target')
    else:
        dir_path = join(dir_path, 'cage_target')

    try:
        os.mkdir(dir_path)
    except:
        pass

    return dir_path


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
    macro_mol.set_position([0,0,0])
    target.set_position([0,0,0])

    # Get the position matrix of the target molecule.
    og_pos_mat = target.position_matrix()

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

        rot_target.set_position_from_matrix(new_pos_mat)

        yield chem.CombineMols(macro_mol.mol, rot_target.mol)

def _c60_rotations(macro_mol, c60, n5fold, n2fold):
    """
    Rotates C60 about its axes of symmetry while placed in `macro_mol`.

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


    macro_mol.set_position([0,0,0])
    c60.set_position([0,0,0])

    # Step 1: Align the 5 membered ring with the z-axis.

    # Find a the ids of atoms in a membered ring.
    g = c60.graph()
    ids = next(x for x in nx.cycle_basis(g) if len(x) == 5)
    # Place the coordinates of those atoms in a matrix.
    ring_matrix = np.matrix([c60.atom_coords(id_) for id_ in ids])

    # Get the centroid of the ring.
    ring_centroid = matrix_centroid(ring_matrix)
    # Align the centroid of the ring with the z-axis.
    c60.set_orientation(ring_centroid, [0,0,1])
    aligned_c60 = copy.deepcopy(c60)

    # Step 2: Get the rotation angles and apply the rotations. Yield
    # the resulting complex.

    # Get the angles of the 5 and 2 fold rotations.
    angles5fold = np.arange(0, 72/180*np.pi, 72/180*np.pi/n5fold)
    angles2fold = np.arange(0, np.pi, np.pi/n2fold)

    for angle5 in angles5fold:
        for angle2 in angles2fold:
            buckyball = copy.deepcopy(aligned_c60)
            buckyball.rotate(angle5, [0,0,1])
            buckyball.rotate(angle2, [0,1,0])
            yield chem.CombineMols(macro_mol.mol, buckyball.mol)
