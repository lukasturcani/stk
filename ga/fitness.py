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
into the `unscaled_fitness` attribute. Next, MMEA copies this value
into the `fitness` attribute at the start of each generation.

So what happens if normalization functions are used?

The first thing to note is that multiple normalization functions can
be applied sequentially. Each normalization function replaces the
previous value in the `fitness` attribute. Normalization functions do
not manipulate or interact with the `unscaled_fitness` attribute in any
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
first argument and this argument should be called ``macro_mol``. It
should also hold a keyword argument called ``logger``. The purpose of
this is to help users identify which arguments are handled
automatically by MMEA and which they need to define in the input file.
The convention is that if the fitness function takes an argument called
``macro_mol`` or ``logger`` they do not have to specify that argument
in the input file.

When defining fitness functions the ``logger`` argument should be
used for logging as a normal logger from the ``logging`` library would.
When running the GA a special logger compatible with multiprocessing
is automatically placed in this argument. It may be useful to define
the logger argument as a keyword argument::

    fit_func(macro_mol, somearg, logger=logging.getLogger(__name__)):
        ...

In this way, if the fitness function is used outside of the GA,
the logger will be provided automatically as well.

A fitness function must return the value which represents the fitness
of the molecule received as an argument. If a fitness function is meant
to be paired with a normalization funtion it can return any value or
object it likes. Just as long as the normalization functions know how
to deal with it and convert it to a number.

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

    macro_mol.progress_params['example_func'] = [mol_energy]

Now a plot showing the change in `mol_energy` across generations will
be made too, along with the plot showing the changes in fitness. In
this case the name of the fitness function was ``example_func``.

What if two things are needed to be kept track of?

    macro_mol.progress_params['example_func'] = [mol_energy, mol_radius]

Great, now a progress plot for each of the variables will be made.

How will the y axes be labelled in each plot?
The decorator `_param_labels()` exists for this.

Let's create a basic outline of a some fitness function:

    @_param_labels('Molecule Energy / J mol-1', 'Mol Radius / m-9')
    def this_is_the_fitness_function(macro_mol, some_param):
        ...
        calculate_stuff()
        ...
        macro_mol.progress_params['this_is_the_fitness_function'] = [
                                                mol_energy, mol_radius]
        ...
        return fitness_value

If this function is used in the GA, a progress plot will be made for
each of the `progress_params` and they will have their y-axes labelled
'Molecule Energy / J mol-1' and 'Molecule Radius / m-9', respectively.

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
from traceback import format_exc

from ..convenience_tools import (matrix_centroid,
                                 FunctionData,
                                 rotation_matrix_arbitrary_axis,
                                 StopLogging, mplogger, FakeLogger)

from ..molecular import (Cage, StructUnit,
                         Energy, optimization, func_key)


logger = logging.getLogger(__name__)


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
    None : NoneType

    """

    # In order for logging to work with multiprocessing properly, each
    # subprocess will log into the que. A thread in the main process
    # will then log.
    m = mp.Manager()
    logq = m.Queue()
    exit_ = StopLogging()
    log_thread = Thread(target=mplogger, args=(logq, logger))
    log_thread.daemon = True
    log_thread.start()

    # Get the fitness function object.
    func = globals()[func_data.name]
    # Make sure it won't raise errors while using multiprocessing.
    p_func = _FitnessFunc(partial(func, **func_data.params), logq)

    # Apply the function to every member of the population, in
    # parallel.
    with mp.get_context('spawn').Pool() as pool:
        evaluated = pool.map(p_func, population)

    # Make sure the cache is updated with the evaluated versions.
    for member in evaluated:
        member.update_cache()

    logq.put((exit_, exit_))
    log_thread.join()


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
    p_func = _FitnessFunc(partial(func, **func_data.params))
    # Apply the function to every member of the population.
    for member in population:
        p_func(member)


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

    def __init__(self, func, logq=None):
        wraps(func)(self)
        self.logq = logq

    def __call__(self, macro_mol, *args,  **kwargs):
        logger = logging.getLogger(__name__)
        if self.logq is not None:
            logger = FakeLogger(self.logq)

        func_name = self.__wrapped__.func.__name__

        # If the fitness function has already been applied to this
        # molecule, return.
        if func_name in macro_mol.unscaled_fitness:
            logger.info('Skipping {}.'.format(macro_mol.name))
            return macro_mol

        try:
            logger.info('Calculating fitness of {}.'.format(
                                                       macro_mol.name))
            val = self.__wrapped__(macro_mol, *args,
                                   **kwargs, logger=logger)

        except Exception as ex:
            val = None
            errormsg = ('Fitness function "{}()" '
                        'failed on molecule "{}".').format(
                        func_name, macro_mol.name)
            logger.error((errormsg+'\n'+format_exc()).strip())

        finally:
            macro_mol.unscaled_fitness[func_name] = val
            return macro_mol


def random_fitness(macro_mol, logger=logger):
    """
    Returns a random fitness value.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule to which a fitness value is to be assigned.

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Returns
    -------
    float
        A random postive number.

    """

    return abs(np.random.normal(50, 20))


@_param_labels('var1', 'var2', 'var3', 'var4')
def random_fitness_vector(macro_mol, logger=logger):
    """
    Returns a size 4 array of random numbers.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule which is to have its fitness calculated.

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Modifies
    --------
    macro_mol.progress_params : dict
        The random numbers are also placed into this attribute.

    Returns
    -------
    numpy.array
        An array holding random numbers.

    """

    # Make a random fitness vector.
    f = abs(np.random.normal(50, 20, 4))
    # This multiplication ensures that the elements of the fitness
    # vector all have different oraders of magnitude and that some
    # are negative.
    f = np.multiply(f, np.array([0.01, 1, 10, -100]))
    macro_mol.progress_params['random_fitness_vector'] = f.tolist()
    return f


def raiser(macro_mol, param1, param2=2, logger=logger):
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

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Returns
    -------
    This function does not return. It only raises.

    Raises
    ------
    Exception
        An exception is always raised.

    """

    raise Exception('Raiser fitness function used.')


@_param_labels('var1', 'var2', 'var3', 'var4')
def partial_raiser(macro_mol, logger=logger):
    """
    Calculates fitness or raises at random.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The molecule having its fitness calculated, maybe.

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Returns
    -------
    numpy.array
        The value of applying random_fitness_vector() to `macro_mol`.

    Raises
    ------
    Exception
        Raised at random.

    """

    if not np.random.choice([0, 1]):
        raise Exception('Partial raiser.')

    r = random_fitness_vector(macro_mol)
    n1 = 'partial_raiser'
    n2 = 'random_fitness_vector'
    macro_mol.progress_params[n1] = macro_mol.progress_params[n2]
    return r


# Provides labels for the progress plotter.
@_param_labels('Cavity Difference ', 'Window Difference ',
               'Asymmetry ', 'Energy per Bond ', 'Precursors Strain',
               'Dihedral Strain')
def cage(macro_mol, dihedral_SMARTS, target_value, pseudoformation_params={
        'func': FunctionData('rdkit', forcefield='mmff')}, logger=logger):
    """
    Returns the fitness vector of a cage.

    The fitness vector consists of the following properties in the
    listed order

        1) `cavity` - the diameter of the cage pore.
        2) `window` - the diameter of the largest cage window.
        3) `asymmetry` - the sum of the size differences of all the
           windows in `macro_mol`.
        4) `eng_per_bond` - The formation energy of `macro_mol` per
           bond made.
        5) `prec_strain` - The mean rmsd between the free building blocks
           and those in the macromolecule.
        6) `dihedral_strain` - The % relative difference between the average
            dihedral angle within the molecule and a target value. The user
            must provide the SMARTS for the dihedral and the target value.

    Notes
    -----
    This function modifies `macro_mol`. It places the calculated
    fitness parameters into :attr:``~.MacroMolecule.progress_params``.

    Parameters
    ----------
    macro_mol : :class:`.Cage`
        The cage whose fitness is to be calculated.

    dihedral_SMARTS : :class:`str`
        The SMARTS code for the dihedral of interest.

    target_value : :class:`float`
        Float representing the target value for the dihedral angle.

    pseudoformation_params : dict (default =
          { 'func' : FunctionData('rdkit', forcefield='uff') })

        This fitness function calculates the formation energy using the
        ``Energy.pseudoformation()`` method. This parameter defines the
        arguments passed to this method via a dictionary. The name of
        the argument is the key and the value of the argument is the
        value.

        Default initialized arguments of Energy.pseudoformation() only
        need to be specified in ``energy_params`` if the user wishes to
        change the default value.

        To see what arguments the ``Energy.pseudoformation()`` method
        requires, try using the  `-h` option:

            python -m mmea -h energy

    logger : :class:`.FakeLogger` or :class:`logging.Logger`, optional
        Used for logging. Not used by this function.


    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holds the fitness vector described in this
        docstring.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """

    # Prevents warnings from getting printed when using
    # multiprocessing.
    warnings.filterwarnings('ignore')

    cavity = macro_mol.cavity_size()
    window = max(macro_mol.windows)
    asymmetry = macro_mol.window_difference()

    logger.debug('Calculating cage energy.')
    e_per_bond = macro_mol.energy.pseudoformation(
                                           **pseudoformation_params)
    e_per_bond /= macro_mol.bonds_made

    prec_strain = macro_mol.bb_distortion()

    dihedral_strain = macro_mol.dihedral_strain(dihedral_SMARTS,
                                                target_value)

    macro_mol.progress_params['cage'] = [cavity, window, asymmetry,
                                         e_per_bond, prec_strain,
                                         dihedral_strain]

    if None in macro_mol.progress_params['cage']:
        raise ValueError(('At least one'
                         ' fitness parameter not calculated.'))

    return np.array([cavity, window, asymmetry, e_per_bond, prec_strain,
                     dihedral_strain])


@_param_labels('Binding Energy', 'Complex Cavity', 'Complex Asymmetry',
               'Complex Strain', 'Cavity', 'Asymmetry', 'Precursors Strain',
               'Dihedral Strain')
def cage_target(macro_mol, target_mol_file, dihedral_SMARTS, target_value,
                efunc, ofunc, rotations=0, logger=logger):
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
        8. strain in the relevant dihedral angles of the cage itself

    Notes
    -----
    This function modifies `macro_mol`. It places the calculated
    fitness parameters into :attr:``~.MacroMolecule.progress_params``.

    Parameters
    ----------
    macro_mol : :class:`.Cage`
        The cage which is to have its fitness calculated.

    target_mol_file : :class:`str`
        The full path of the ``.mol`` file hodling the target molecule
        placed inside the cage.

    dihedral_SMARTS : :class:`str`
        The SMARTS code for the dihedral of interest.

    target_value : :class:`float`
        Float representing the target value for the dihedral angle.

    efunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the energy
        function used to calculate energies.

    ofunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the optimization
        function to be run on the generated complexes.

    rotations : :class:`int`, optional
        The number of times the target should be randomly rotated
        within the cage cavity in order to find the most stable
        conformation.

    logger : :class:`.FakeLogger` or :class:`logging.Logger`, optional
        Used for logging. Not used by this function.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holds the fitness vector described in this
        docstring.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """

    return _cage_target('cage_target', macro_mol,
                        target_mol_file, dihedral_SMARTS,
                        target_value, efunc, ofunc,
                        FunctionData('_generate_complexes',
                                     number=rotations+1),
                        logger)


@_param_labels('Binding Energy', 'Complex Cavity', 'Complex Asymmetry',
               'Complex Strain', 'Cavity', 'Asymmetry', 'Precursors Strain',
               'Dihedral Strain')
def cage_c60(macro_mol, target_mol_file, dihedral_SMARTS, target_value,
             efunc, ofunc, n5fold, n2fold, logger=logger):
    """
    Calculates the fitness vector of a cage / C60 complex.

    The difference between this function and :func:`cage_target` is
    that the rotations are specifically aimed at sampling C60 entirely
    and systematically. Rather than the random sampling of the other
    function.

    The function returns a fitness vector consisting of:

        1. binding energy
        2. cavity of cage in complex
        3. asymmetry of cage in complex
        4. strain of cage in complex
        5. cavity of cage by itself
        6. asymmetry of cage by itself
        7. strain of cage by itself
        8. strain in the relevant dihedral angles of the cage itself

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

    dihedral_SMARTS : :class:`str`
        The SMARTS code for the dihedral of interest.

    target_value : :class:`float`
        Float representing the target value for the dihedral angle.

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

    logger : :class:`.FakeLogger` or :class:`logging.Logger`, optional
        Used for logging. Not used by this function.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holds the fitness vector described in this
        docstring.

    Raises
    ------
    :class:`ValueError`
        If the calculation of a fitness parameter fails.

    """
    return _cage_target('cage_c60', macro_mol,
                        target_mol_file, efunc, ofunc,
                        FunctionData('_c60_rotations',
                                     n5fold=n5fold,
                                     n2fold=n2fold),
                        logger)


def _cage_target(func_name, macro_mol, target_mol_file, dihedral_SMARTS,
                 target_value, efunc, ofunc, rotation_func, logger):
    """
    A general fitness function for calculating fitness of complexes.

    This function should be inherited by other fitness functions which
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
        8. strain in the relevant dihedral angles of the cage itself

    Notes
    -----
    This function modifies `macro_mol`. It places the calculated
    fitness parameters into :attr:`~.MacroMolecule.progress_params`.

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

    dihedral_SMARTS : :class:`str`
        The SMARTS code for the dihedral angle of interest.

    target_value : :class:`float`
        Float representing the target value for the dihedral angle.

    efunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the energy
        function used to calculate energies.

    ofunc : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the optimization
        function to be run on the generated complexes.

    rotation_func : :class:`.FunctionData`
        A :class:`.FunctionData` object representing the rotation
        function to be used.

    logger : :class:`.FakeLogger` or :class:`logging.Logger`
        Used for logging. Not used by this function.

    Returns
    -------
    :class:`numpy.ndarray`
        The numpy array holds the fitness vector described in this
        docstring.

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
                                                target_value)
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
                     cmplx_cavity, cmplx_asymmetry, cmplx_strain,
                     cavity, asymmetry, prec_strain, dihedral_strain])


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
