"""
Module for defining fitness functions.

How fitness values are calculated.
----------------------------------

The basic outline of the fitness calculation procedure for a single
molecule is as follows,

.. code-block:: python

    # 1. Start with a MacroMolecule instance.
    mol = Polymer([monomer1, monomer2], SomeTopology())

    # This new instance has had no fitness functions applied to it
    # and therefore its "fitness" and "unscaled_fitness" attributes
    # are empty.
    mol.fitness  # None
    mol.unscaled_fitness  # {}

    # 2. Apply a fitness function to the MacroMolecule.
    mol.unscaled_fitness['random_fitness'] = random_fitness(mol)

    # This adds the unscaled fitness value into the "unscaled_fitness"
    # attribute, but does not change the "fitness" attribute.
    mol.fitness  # None
    mol.unscaled_fitness  # {'random_fitness': 34}

    # 3.Copy the value from "unscaled_fitness" into "fitness".
    mol.fitness = deepcopy(mol.unscaled_fitness['random_fitness'])

    mol.fitness  # 34
    mol.unscaled_fitness  # {'random_fitness': 34}

    # 4. Apply any number of normalization functions to the fitness
    #    value. Normalization functions can be skipped altogether.

    # Normalization functions act on a population of MacroMolecules,
    # i.e. they normalize the fitness values across the population. So
    # first step is to put the molecule into a population.
    pop = GAPopulation(mol)

    # Apply the normalization function.
    norm_func(pop)

    mol.fitness  # 12
    mol.unscaled_fitness  # {'random_fitness': 34}

    # Multiple normalization functions may be applied.
    norm_func2(pop)

    mol.fitness  # 200
    mol.unscaled_fitness  # {'random_fitness': 34}

When running the GA, this process is much simplified,

.. code-block:: python

    # Start with a population of molecules, representing a GA
    # generation.
    pop = GAPopulation(mol1, mol2, mol3)

    mol1.fitness  # None
    mol1.unscaled_fitness  # {}

    mol2.fitness  # None
    mol2.unscaled_fitness  # {}

    mol3.fitness  # None
    mol3.unscaled_fitness  # {}

    # Calculate the unscaled fitness values of all molecules.
    pop.calculate_member_fitness()

    mol1.fitness  # None
    mol1.unscaled_fitness  # {'random_fitness': 12}

    mol2.fitness  # None
    mol2.unscaled_fitness  # {'random_fitness': 34}

    mol3.fitness  # None
    mol3.unscaled_fitness  # {'random_fitness': 22}

    # Apply the normalization functions.
    pop.normalize_fitness_values()

    mol1.fitness  # 3
    mol1.unscaled_fitness  # {'random_fitness': 12}

    mol2.fitness  # 1
    mol2.unscaled_fitness  # {'random_fitness': 34}

    mol3.fitness  # 2
    mol3.unscaled_fitness  # {'random_fitness': 22}

The method :meth:`.GAPopulation.normalize_fitness_values` automatically
copies the fitness value from :attr:`.MacroMolecule.unscaled_fitness`
into :attr:`.MacroMolecule.fitness` before applying any normalization
functions. A fitness function is only applied once per molecule,
because the result depends only on the molecule. However, normalization
functions are re-applied every generation. This is because for some
normalization functions, the returned value for a molecule may also
depend on the other molecules in the population. At each generation,
a fresh copy of the value from :attr:`.MacroMolecule.unscaled_fitness`
is made into :attr:`.MacroMolecule.fitness`. This means that the
normalization functions always start from the unscaled fitness value.

Extending mtk: Adding fitness functions.
----------------------------------------

To add a new fitness function simply write it as a function in this
module. It will need to take the :class:`.MacroMolecule` instance as
its first argument and this argument should be called `macro_mol`. The
purpose of this is to help users identify which arguments are handled
automatically by ``mtk`` and which they need to define in the input
file. The convention is that if the fitness function takes an argument
called `macro_mol`, they do not have to specify that argument in the
input file.

A fitness function must return the value which represents the fitness
of the molecule `macro_mol`. If a fitness function is meant to be
paired with a normalization function it can return any value or
object it likes. Just as long as the normalization functions know how
to deal with it and convert it to a number.

.. code-block:: python

    def random_fitness(macro_mol):
        return abs(np.random.normal(50, 20))

Obviously, this is a just toy example but all the key components of
fitness functions are there. More complex fitness functions can take
an arbitrary number arguments, and will likely use molecular properties
to calculate a fitness value. Here is a fitness function which rewards
molecules for having a desired size.

.. code-block:: python

    def size(macro_mol, desired_size):
        return abs(macro_mol.max_diameter() - desired_size)

A fitness function may be complex and may not fit neatly into a single
function. For example, the :func:`cage_target` needs to call
:func:`_generate_complexes` in order to sample various conformations
before outputting a fitness value. This is fine. Define helper
functions such as :func:`_generate_complexes` within this module but
make sure they are private. This means that names of helper functions
begin with a leading underscore.

.. _plotting-note:

A note on plotting.
-------------------

As mentioned before, some fitness functions may be complex and as a
result manipulate all sorts of data. Typically, in order to measure the
progress of a GA, the fitness values in the population are tracked
across generations. However, let's say that some hypothetical fitness
function also calculates the energies of molecules. It may be quite
interesting plot the evolution of energies across generations too, even
if the energy is not directly reflect in the final fitness value. If
this is the case, the fitness function may assign to the
:attr:`.MacroMolecule.progress_params` attribute of `macro_mol`,

.. code-block:: python

    def example_func(macro_mol):
        mol_energy = macro_mol.energy.some_energy_func()
        macro_mol.progress_params['example_func'] = [mol_energy]
        ...

Now a plot showing the change in ``mol_energy`` across generations will
be made too, along with the plot showing the changes in fitness.

What if two things are needed to be kept track of? Simple,

.. code-block:: python

    def example_func(macro_mol):
        mol_energy = macro_mol.energy.some_energy_func()
        mol_radius = macro_mol.max_diamater() / 2
        macro_mol.progress_params['example_func'] = [mol_energy,
                                                     mol_radius]
        ...

Great, now a separate progress plot ``mol_energy`` and ``mol_radius``
will be made.

How will the y-axes be labelled in each plot? The decorator
:func:`_param_labels` exists for this.

Let's create a basic outline of a some fitness function:

.. code-block:: python

    @_param_labels('Molecule Energy / J mol-1', 'Mol Radius / m-9')
    def some_fitness_fn(macro_mol, some_param):
        ...
        calculate_stuff()
        ...
        macro_mol.progress_params['some_fitness_fn'] = [mol_energy,
                                                        mol_radius]
        ...
        return fitness_value

If this function is used in the GA, a progress plot will be made for
each variable placed in :attr:`.MacroMolecule.progress_params` and
they will have their y-axes labelled ``'Molecule Energy / J mol-1'``
and ``'Molecule Radius / m-9'``, respectively.

"""

import numpy as np
import rdkit.Chem.AllChem as rdkit
import copy
import os
import warnings
from functools import partial, wraps
import networkx as nx
import multiprocessing as mp
from collections import Counter
from os.path import join
from uuid import uuid4
import logging
from threading import Thread

from ..convenience_tools import (matrix_centroid,
                                 FunctionData,
                                 rotation_matrix_arbitrary_axis,
                                 daemon_logger,
                                 logged_call)

from ..molecular import Cage, StructUnit, Energy, func_key
from .. import optimization

logger = logging.getLogger(__name__)


def _calc_fitness(func_data, population, processes):
    """
    Calculates the fitness values of all members of a population.

    Parameters
    ----------
    func_data : :class:`.FunctionData`
        A :class:`.FunctionData` instance representing the chosen
        fitness function and any additional parameters it may require.

    population : :class:`.GAPopulation`
        The population whose members must have their fitness
        calculated.

    processes : :class:`int`
        The number of parallel processes to create.

    Returns
    -------
    None : :class:`NoneType`

    """

    manager = mp.Manager()
    logq = manager.Queue()
    log_thread = Thread(target=daemon_logger, args=(logq, ))
    log_thread.start()

    # Get the fitness function object.
    func = globals()[func_data.name]
    # Make sure it won't raise errors while using multiprocessing.
    p_func = _FitnessFunc(partial(func, **func_data.params))

    # Apply the function to every member of the population, in
    # parallel.
    with mp.get_context('spawn').Pool(processes) as pool:
        evaluated = pool.starmap(logged_call,
                                 ((logq, p_func, mem) for
                                  mem in population))

    # Make sure the cache is updated with the evaluated versions.
    for member in evaluated:
        member.update_cache()

    logq.put(None)
    log_thread.join()


def _calc_fitness_serial(func_data, population):
    """
    Calculates the fitness values of all members of a population.

    Parameters
    ----------
    func_data : :class:`.FunctionData`
        A :class:`.FunctionData` instance representing the chosen
        fitness function and any additional parameters it may require.

    population : :class:`.GAPopulation`
        The population whose members must have their fitness
        calculated.

    Returns
    -------
    None : :class:`NoneType`

    """

    # Get the fitness function object.
    func = globals()[func_data.name]
    p_func = _FitnessFunc(partial(func, **func_data.params))
    # Apply the function to every member of the population.
    for member in population:
        p_func(member)


def _param_labels(*labels):
    """
    Adds the :attr:`param_labels` attribute to a fitness function.

    The point of this decorator is described in :ref:`plotting-note`.

    Parameters
    ----------
    *labels : :class:`str`
        This function takes an arbitrary number of strings. The strings
        are the y-axis labels for each graph made per variable in
        :attr:`.MacroMolecule.progress_params`. The order of strings
        must correspond to the order of variables placed into
        :attr:`.MacroMolecule.progress_params`.

    Returns
    -------
    :class:`function`
        The decorated function, with the attribute :attr:`param_labels`
        added. :attr:`param_labels` holds the strings provided in
        `*labels`.

    Examples
    --------

    .. code-block:: python

        @_param_labels('ONE', 'TWO', 'THREE')
        def some_func(...):
            return 1

        some_func.param_labels  # ['ONE', 'TWO', 'THREE']

    """

    def add_labels(func):
        func.param_labels = labels
        return func

    return add_labels


class _FitnessFunc:
    """
    A decorator for fitness functions.

    This decorator is applied to all fitness functions automatically in
    :func:`_calc_fitness`. It should not be applied explicitly when
    defining the functions.

    The decorator prevents fitness functions from raising if
    they fail (necessary for ``multiprocessing`` compatibility),
    prevents them from being run twice on the same molecule and stores
    the value returned by them in
    :attr:`.MacroMolecule.unscaled_fitness`.

    """

    def __init__(self, func):
        """
        Initializes a :class:`_FitnessFunc` instance.

        Parameters
        ----------
        func : :class:`function`
            The fitness function to be decorated.

        """

        wraps(func)(self)

    def __call__(self, macro_mol):
        """
        Decorates and calls the fitness function.

        Parameters
        ----------
        macro_mol : :class:`.MacroMolecule`
            The molecule to have its fitness calculated.

        Returns
        -------
        :class:`.MacroMolecule`
            `macro_mol` with its fitness calculated.

        """

        func_name = self.__wrapped__.func.__name__

        # If the fitness function has already been applied to this
        # molecule, return.
        if func_name in macro_mol.unscaled_fitness:
            logger.info(f'Skipping {macro_mol.name}.')
            return macro_mol

        try:
            logger.info(f'Calculating fitness of {macro_mol.name}.')
            val = self.__wrapped__(macro_mol)

        except Exception as ex:
            val = None
            errormsg = (f'Fitness function "{func_name}()" '
                        f'failed on molecule "{macro_mol.name}".')
            logger.error(errormsg, exc_info=True)

        finally:
            macro_mol.unscaled_fitness[func_name] = val
            return macro_mol


def random_fitness(macro_mol):
    """
    Returns a random fitness value.

    Parameters
    ----------
    macro_mol : :class:`.MacroMolecule`
        The molecule for which a fitness value is to be calculated.

    Returns
    -------
    :class:`float`
        A random postive number.

    """

    return abs(np.random.normal(50, 20))


@_param_labels('var1', 'var2', 'var3', 'var4')
def random_fitness_vector(macro_mol):
    """
    Returns a 4 element array of random numbers.

    Notes
    -----
    This function places the array into
    :attr:`~.MacroMolecule.progress_params` of `macro_mol`.

    Parameters
    ----------
    macro_mol : :class:`.MacroMolecule`
        The molecule for which a fitness value is to be calculated.

    Returns
    -------
    :class:`numpy.ndarray`
        An array holding the 4 random numbers.

    """

    # Make a random fitness vector.
    f = abs(np.random.normal(50, 20, 4))
    # This multiplication ensures that the elements of the fitness
    # vector all have different oraders of magnitude and that some
    # are negative.
    f = np.multiply(f, np.array([0.01, 1, 10, -100]))
    macro_mol.progress_params['random_fitness_vector'] = f.to_list()
    return f


def raiser(macro_mol, param1, param2=2):
    """
    Doens't calculate a fitness value, raises an error instead.

    This function is used for tests to ensure that when fitness
    functions raise errors, they are handeled correctly.

    Parameters
    ---------
    param1 : :class:`object`
        Dummy parameter, does nothing.

    param2 : :class:`object`, optional
        Dummy keyword parameter, does nothing.

    Returns
    -------
    None : :class:`NoneType`
        This function does not return. It only raises.

    Raises
    ------
    :class:`Exception`
        An exception is always raised.

    """

    raise Exception('Raiser fitness function used.')


@_param_labels('var1', 'var2', 'var3', 'var4')
def partial_raiser(macro_mol):
    """
    Calculates fitness or raises at random.

    Parameters
    ----------
    macro_mol : :class:`.MacroMolecule`
        The molecule having its fitness calculated.

    Returns
    -------
    :class:`numpy.ndarray`
        The result of applying :func:`random_fitness_vector` to
        `macro_mol`.

    Raises
    ------
    :class:`Exception`
        Raised at random.

    """

    if np.random.choice([0, 1]):
        raise Exception('Partial raiser.')

    r = random_fitness_vector(macro_mol)
    n1 = 'partial_raiser'
    n2 = 'random_fitness_vector'
    macro_mol.progress_params[n1] = macro_mol.progress_params[n2]
    return r


# Provides labels for the progress plotter.
@_param_labels('Cavity Difference',
               'Window Difference',
               'Asymmetry',
               'Energy per Bond',
               'Precursors Strain',
               'Dihedral Strain')
def cage(macro_mol,
         pseudoformation_params={'func': FunctionData('rdkit',
                                                      forcefield='mmff')},
         dihedral_SMARTS='',
         target_dihedral=180):
    """
    Returns the fitness vector of a cage.

    The fitness vector consists of the following properties in the
    listed order

        1. `cavity` - the diameter of the cage pore.
        2. `window` - the diameter of the largest cage window.
        3. `asymmetry` - the sum of the size differences of all the
           windows in `macro_mol`.
        4. `eng_per_bond` - The formation energy of `macro_mol` per
           bond made.
        5. `prec_strain` - The mean rmsd between the free building
           block and those in the macromolecule.
        6. `dihedral_strain` - The % relative difference between the
           average dihedral angle within the molecule and a target
           value. The user must provide the SMARTS for the dihedral and
           the target value.

    Parameters
    ----------
    macro_mol : :class:`.Cage`
        The cage whose fitness is to be calculated.

    dihedral_SMARTS : :class:`str`, optional
        The SMARTS code for the dihedral of interest.

    target_dihedral : :class:`float`, optional
        The target value for the dihedral angle.

    pseudoformation_params : dict, optional
        This fitness function calculates the formation energy using
        :meth:`.Energy.pseudoformation`. This parameter defines the
        arguments passed to this method via a dictionary. The name of
        the argument is the key and the value of the argument is the
        value.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holding the fitness vector.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """

    # Prevents warnings from getting printed when using
    # multiprocessing.
    warnings.filterwarnings('ignore')

    cavity = macro_mol.cavity_size()
    window = max(macro_mol.windows())
    asymmetry = macro_mol.window_difference()

    logger.debug('Calculating cage energy.')
    e_per_bond = macro_mol.energy.pseudoformation(
                                           **pseudoformation_params)
    e_per_bond /= macro_mol.bonds_made

    prec_strain = macro_mol.bb_distortion()

    dihedral_strain = macro_mol.dihedral_strain(dihedral_SMARTS,
                                                target_dihedral)

    macro_mol.progress_params['cage'] = [cavity,
                                         window,
                                         asymmetry,
                                         e_per_bond,
                                         prec_strain,
                                         dihedral_strain]

    if None in macro_mol.progress_params['cage']:
        raise ValueError(('At least one'
                         ' fitness parameter not calculated.'))

    return np.array([cavity,
                     window,
                     asymmetry,
                     e_per_bond,
                     prec_strain,
                     dihedral_strain])


@_param_labels('Binding Energy',
               'Complex Cavity',
               'Complex Asymmetry',
               'Complex Strain',
               'Cavity',
               'Asymmetry',
               'Precursors Strain',
               'Dihedral Strain')
def cage_target(macro_mol,
                target_mol_file,
                efunc,
                ofunc,
                dihedral_SMARTS="",
                target_value=180,
                rotations=0):
    """
    Returns the fitness vector of a cage / target complex.

    The target is randomly rotated inside the cage's cavity and the
    most stable conformation found is used.

    The function returns a fitness vector consisting of:

        1. binding energy
        2. cavity of cage in complex
        3. asymmetry of cage in complex
        4. strain of cage in complex
        5. cavity of cage by itself
        6. asymmetry of cage by itself
        7. strain of cage by itself
        8. strain in select dihedral angles of the cage by itself

    Notes
    -----
    This function modifies `macro_mol`. It places the calculated
    fitness parameters into :attr:`~.MacroMolecule.progress_params`.

    Parameters
    ----------
    macro_mol : :class:`.Cage`
        The cage which is to have its fitness calculated.

    target_mol_file : :class:`str`
        The full path of the ``.mol`` file hodling the target molecule
        placed inside the cage.

    efunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the energy
        function used to calculate energies.

    ofunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the optimization
        function to be run on the generated complexes.

    dihedral_SMARTS : :class:`str`, optional
        The SMARTS code for the dihedral of interest.

    target_value : :class:`float`, optional
        A number representing the target value for the dihedral angle.

    rotations : :class:`int`, optional
        The number of times the target should be randomly rotated
        within the cage cavity in order to find the most stable
        conformation.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holding the fitness vector.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """

    return _cage_target('cage_target',
                        macro_mol,
                        target_mol_file,
                        efunc,
                        ofunc,
                        FunctionData('_generate_complexes',
                                     number=rotations+1),
                        dihedral_SMARTS,
                        target_value)


@_param_labels('Binding Energy',
               'Complex Cavity',
               'Complex Asymmetry',
               'Complex Strain',
               'Cavity',
               'Asymmetry',
               'Precursors Strain',
               'Dihedral Strain')
def cage_c60(macro_mol,
             target_mol_file,
             efunc,
             ofunc,
             n5fold,
             n2fold,
             dihedral_SMARTS='',
             target_dihedral=180):
    """
    Calculates the fitness vector of a cage / C60 complex.

    The difference between this function and :func:`cage_target` is
    that the rotations are specifically aimed at sampling C60 entirely
    and systematically. Rather than the random sampling used by the
    other function.

    The function returns a fitness vector consisting of:

        1. binding energy
        2. cavity of cage in complex
        3. asymmetry of cage in complex
        4. strain of cage in complex
        5. cavity of cage by itself
        6. asymmetry of cage by itself
        7. strain of cage by itself
        8. strain in select dihedral angles of the cage by itself

    Parameters
    ----------
    macro_mol : :class:`.Cage`
        The cage which is to have its fitness calculated.

    target_mol_file : :class:`str`
        The full path of the ``.mol`` file hodling the target molecule
        placed inside the cage.

    efunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the energy
        function used to calculate energies.

    ofunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the optimization
        function to be run on the generated complexes.

    n5fold : :class:`int`
        The number of rotations along the 5-fold axis of symmetry.

    n2fold : :class:`int`
        The number of rotations along the 2 fold axis of symmetry per
        rotation along the 5-fold axis.

    dihedral_SMARTS : :class:`str`, optional
        The SMARTS code for the dihedral of interest.

    target_dihedral : :class:`float`, optional
        The target value for the dihedral angle.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holding the fitness vector.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """
    return _cage_target('cage_c60',
                        macro_mol,
                        target_mol_file,
                        efunc,
                        ofunc,
                        FunctionData('_c60_rotations',
                                     n5fold=n5fold,
                                     n2fold=n2fold),
                        dihedral_SMARTS,
                        target_dihedral)


def _cage_target(func_name,
                 macro_mol,
                 target_mol_file,
                 efunc,
                 ofunc,
                 rotation_func,
                 dihedral_SMARTS,
                 target_dihedral):
    """
    A general fitness function for calculating fitness of complexes.

    This function should be wrapped by other fitness functions which
    define their own rotation function. For example :func:`cage_c60`
    and :func:`cage_target`.

    The function returns a fitness vector consisting of:

        1. binding energy
        2. cavity of cage in complex
        3. asymmetry of cage in complex
        4. strain of cage in complex
        5. cavity of cage by itself
        6. asymmetry of cage by itself
        7. strain of cage by itself
        8. strain in select dihedral angles of the cage by itself

    Parameters
    ----------
    func_name : :class:`str`
        The name of the external fitness function calling this one.
        Used for the key in :attr:`~.MacroMolecule.progress_params`.

    macro_mol : :class:`.Cage`
        The cage which is to have its fitness calculated.

    target_mol_file : :class:`str`
        The full path of the ``.mol`` file hodling the target molecule
        placed inside the cage.

    efunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the energy
        function used to calculate energies.

    ofunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the optimization
        function to be run on the generated complexes.

    rotation_func : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the rotation
        function to be used.

    dihedral_SMARTS : :class:`str`, optional
        The SMARTS code for the dihedral of interest.

    target_dihedral : :class:`float`, optional
        The target value for the dihedral angle.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holding the fitness vector.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """

    warnings.filterwarnings('ignore')

    # Create a folder to hold all the complexes generated by this
    # function.
    folder_path = _make_cage_target_folder()
    # Transform the FunctionData instances into functions.
    optfunc = getattr(optimization, ofunc.name)
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
    rdkit_complexes = rot_func(unopt_macro_mol, target,
                               **rotation_func.params)

    # Optimize the strcuture of the cage/target complexes.
    macromol_complexes = []
    logger.debug('Optimizing complex structures.')
    for i, complex_ in enumerate(rdkit_complexes):
        # In order to use the optimization functions, first the data
        # is loaded into a ``Cage`` instance and its .mol
        # file is written to the disk.
        mm_complex = Cage.__new__(Cage)
        mm_complex.mol = complex_
        mm_complex.name = name + '_COMPLEX_{0}'.format(i)
        mm_complex.optimized = False
        mm_complex.energy = Energy(mm_complex)
        mm_complex.topology = macro_mol.topology
        mm_complex.building_blocks = macro_mol.building_blocks
        mm_complex.bonder_ids = macro_mol.bonder_ids
        optfunc(mm_complex, **ofunc.params)
        macromol_complexes.append(mm_complex)

    # Calculate the energy of the complex and compare to the
    # individual energies. If more than complex was made, use the
    # most stable version.
    mm_energy = getattr(macro_mol.energy,
                        efunc.name)(**efunc.params)
    target_energy = getattr(target.energy, efunc.name)(**efunc.params)

    energy_separate = mm_energy + target_energy

    logger.debug('Calculating complex energies.')
    min_eng_cmplx = min(macromol_complexes,
                        key=lambda x:
                        getattr(x.energy, efunc.name)(**efunc.params))

    # Write the most stable complex to a file.
    min_eng_cmplx.write(join(folder_path, min_eng_cmplx.name+'.mol'))

    ekey = func_key(getattr(Energy, efunc.name),
                    (min_eng_cmplx.energy, ), efunc.params)

    binding_energy = (min_eng_cmplx.energy.values[ekey] -
                      energy_separate)

    frag1, frag2 = rdkit.GetMolFrags(min_eng_cmplx.mol,
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

    cmplx_cage = Cage.__new__(Cage)
    cmplx_cage.building_blocks = list(macro_mol.building_blocks)
    cmplx_cage.bonder_ids = list(macro_mol.bonder_ids)
    cmplx_cage.fragments = macro_mol.fragments
    cmplx_cage.fg_ids = set(macro_mol.fg_ids)
    cmplx_cage.mol = cmplx_cage_mol
    cmplx_cage.topology = macro_mol.topology
    cmplx_cage.name = min_eng_cmplx.name + '_no_target'
    # Calculate fitness parameters of cage in complex.
    cmplx_cavity = cmplx_cage.cavity_size()
    cmplx_asymmetry = cmplx_cage.window_difference()
    cmplx_strain = cmplx_cage.bb_distortion()

    # Write the cage without the target to a file.
    cmplx_cage.write(join(folder_path, cmplx_cage.name+'.mol'))

    # Calculate fitness parameters of cage by itself.
    cavity = macro_mol.cavity_size()
    asymmetry = macro_mol.window_difference()
    prec_strain = macro_mol.bb_distortion()
    dihedral_strain = macro_mol.dihedral_strain(dihedral_SMARTS,
                                                target_dihedral)
    macro_mol.progress_params[func_name] = [binding_energy,
                                            cmplx_cavity,
                                            cmplx_asymmetry,
                                            cmplx_strain,
                                            cavity,
                                            asymmetry,
                                            prec_strain,
                                            dihedral_strain]

    if None in macro_mol.progress_params[func_name]:
        raise ValueError(('At least one'
                         ' fitness parameter not calculated.'),
                         macro_mol.progress_params[func_name])

    return np.array([binding_energy,
                     cmplx_cavity,
                     cmplx_asymmetry,
                     cmplx_strain,
                     cavity,
                     asymmetry,
                     prec_strain,
                     dihedral_strain])


def _make_cage_target_folder():
    """
    Creates a folder to store molecules made by :func:`_cage_target`.

    The function creates a folder called ``cage_target``.
    Inside will be any complexes formed by the :func:`_cage_target`.
    The folder will be placed in the current working directory, or
    if the GA is running, 1 above the current working dirctory. This
    prevents the generated molecules from being cleaned up by the GA
    when it's finished.

    Returns
    -------
    :class:`str`
        The path of the ``cage_target`` folder.

    """

    dir_path = os.getcwd()
    if join('output', 'scratch') in dir_path:
        dir_path = dir_path.replace('scratch', 'cage_target')
    else:
        dir_path = join(dir_path, 'cage_target')

    try:
        os.mkdir(dir_path)
    except Exception:
        pass

    return dir_path


def _generate_complexes(macro_mol, target, number=1):
    """
    Yields ``rdkit`` molecules of cage / target complexes.

    If multiple complexes are returned, they will be different via a
    random rotation accross the x, y and z axes.

    Parameters
    ----------
    macro_mol : :class:`.Cage`
        The cage used to form the complex.

    target : :class:`.StructUnit`
        The target used to form the complex.

    number : :class:`int`, optional
        The number of complexes to be returned.

    Yields
    ------
    :class:`rdkit.Chem.rdchem.Mol`
        An ``rdkit`` instance holding the cage / target complex.

    """

    # First place both the target and cage at the origin.
    macro_mol.set_position([0, 0, 0])
    target.set_position([0, 0, 0])

    # Get the position matrix of the target molecule.
    og_pos_mat = target.position_matrix()

    # Carry out every rotation and yield a complex for each case.
    for i in range(number):
        rot_target = copy.deepcopy(target)

        rot1 = np.random.rand() * 2*np.pi
        rot2 = np.random.rand() * 2*np.pi
        rot3 = np.random.rand() * 2*np.pi

        rot_mat1 = rotation_matrix_arbitrary_axis(rot1, [1, 0, 0])
        rot_mat2 = rotation_matrix_arbitrary_axis(rot2, [0, 1, 0])
        rot_mat3 = rotation_matrix_arbitrary_axis(rot3, [0, 0, 1])

        new_pos_mat = np.dot(rot_mat1, og_pos_mat)
        new_pos_mat = np.dot(rot_mat2, new_pos_mat)
        new_pos_mat = np.dot(rot_mat3, new_pos_mat)

        rot_target.set_position_from_matrix(new_pos_mat)

        yield rdkit.CombineMols(macro_mol.mol, rot_target.mol)


def _c60_rotations(macro_mol, c60, n5fold, n2fold):
    """
    Rotates C60 about its axes of symmetry while placed in `macro_mol`.

    Parameters
    ----------
    macro_mol : :class:`.MacroMolecule`
        The cage which should have C60 placed inside it.

    c60 : :class:`.StructUnit`
        A StructUnit instance of C60.

    n5fold : :class:`int`
        The number of rotations along the 5-fold axis of symmetry.

    n2fold : :class:`int`
        The number of rotations along the 2-fold axis of symmetry per
        rotation along the 5-fold axis.

    Yields
    ------
    :class:`rdkit.Chem.rdchem.Mol`
        An ``rdkit`` instance holding the cage / C60 complex.

    """

    macro_mol.set_position([0, 0, 0])
    c60.set_position([0, 0, 0])

    # Step 1: Align the 5 membered ring with the z-axis.

    # Find a the ids of atoms in a membered ring.
    g = c60.graph()
    ids = next(x for x in nx.cycle_basis(g) if len(x) == 5)
    # Place the coordinates of those atoms in a matrix.
    ring_matrix = np.matrix([c60.atom_coords(id_) for id_ in ids])

    # Get the centroid of the ring.
    ring_centroid = matrix_centroid(ring_matrix)
    # Align the centroid of the ring with the z-axis.
    c60.set_orientation(ring_centroid, [0, 0, 1])
    aligned_c60 = copy.deepcopy(c60)

    # Step 2: Get the rotation angles and apply the rotations. Yield
    # the resulting complex.

    # Get the angles of the 5 and 2 fold rotations.
    angles5fold = np.arange(0, 72/180*np.pi, 72/180*np.pi/n5fold)
    angles2fold = np.arange(0, np.pi, np.pi/n2fold)

    for angle5 in angles5fold:
        for angle2 in angles2fold:
            buckyball = copy.deepcopy(aligned_c60)
            buckyball.rotate(angle5, [0, 0, 1])
            buckyball.rotate(angle2, [0, 1, 0])
            yield rdkit.CombineMols(macro_mol.mol, buckyball.mol)
