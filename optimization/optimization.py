"""
Defines optimization functions.

Extending MMEA: Adding optimization functions
---------------------------------------------
New optimization functions are added by writing them into this module.
The only requirement is that the first argument is ``mol`` and
there is a keyword argument ``logger``. This allows users to identify
which arguments are handled automatically by MMEA and which need to be
defined in the input file. The convention is that if the optimization
function takes the arguments ``mol`` and ``logger``, the user
does not have to specify them argument in the input file.

When defining optimization functions the ``logger`` argument should be
used for logging as a normal logger from the ``logging`` library would.
When running the GA a special logger compatible with multiprocessing
is automatically placed in this argument. It may be useful to define
the logger argument as a keyword argument::

    opt_func(mol, somearg, logger=logging.getLogger(__name__)):
        ...

In this way, if the optimization function is used outside of the GA,
the logger will be provided automatically as well.

An optimization function should update the rdkit molecule in the `mol`
attribute. The return values of optimization functions are discarded.

Optimizations can be complicated. If the use of helper functions is
required make sure that they are private, ie that their names begin
with a leading underscore. In the event that the optimization is so
complex that it requires its own module or file, place it in the same
folder as this file. Then import the optimization function into this
file. See ``macromodel.py`` as an example. Make sure that only the
optimization functions are imported back into this file, not any of the
helper functions or classes.

"""

import rdkit.Chem.AllChem as rdkit
import multiprocessing as mp
from functools import partial, wraps
import numpy as np
import logging
from threading import Thread

from .macromodel import macromodel_opt, macromodel_cage_opt

from ..convenience_tools import daemon_logger, logged_call


logger = logging.getLogger(__name__)


def _optimize_all(func_data, population, processes):
    """
    Run opt function on all population members in parallel.

    Parameters
    ----------
    func_data : FunctionData
        The ``FunctionData`` object which represents the chosen
        optimization function. This function should be defined within
        this module. The ``FunctionData`` object also holds any
        additional parameters the optimization function may need.

    population : Population
        The ``Population`` instance who's members must be optimized.

    processes : :class:`int`
        The number of parallel processes to create.

    Returns
    -------
    None : NoneType

    """

    manager = mp.Manager()
    logq = manager.Queue()
    log_thread = Thread(target=daemon_logger, args=(logq, ))
    log_thread.start()

    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the
    # module.
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may
    # require.
    p_func = _OptimizationFunc(partial(func, **func_data.params))

    # Apply the function to every member of the population, in
    # parallel.
    with mp.get_context('spawn').Pool(processes) as pool:
        optimized = pool.starmap(logged_call,
                                 ((logq, p_func, mem) for
                                  mem in population))
    # Make sure the cache is updated with the optimized versions.
    for member in optimized:
        member.update_cache()

    logq.put(None)
    log_thread.join()


def _optimize_all_serial(func_data, population):
    """
    Run opt function on all population members sequentially.

    Parameters
    ----------
    func_data : FunctionData
        The ``FunctionData`` object which represents the chosen
        optimization function. This function should be defined within
        this module. The ``FunctionData`` object also holds any
        additional parameters the optimization function may need.

    population : Population
        The ``Population`` instance who's members must be optimized.

    Returns
    -------
    None : NoneType


    """

    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the
    # module.
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may
    # require.
    p_func = _OptimizationFunc(partial(func, **func_data.params))

    # Apply the function to every member of the population.
    for member in population:
        p_func(member)


class _OptimizationFunc:
    """
    A decorator for optimziation functions.

    This decorator is applied to all optimization functions
    automatically in _optimize_all(). It should not be applied
    explicitly when defining the functions.

    This decorator prevents optimization functions from raising if
    they fail (necessary for multiprocessing) and prevents them from
    being run twice on the same molecule.

    """

    def __init__(self, func):
        wraps(func)(self)

    def __call__(self, mol):

        if mol.optimized:
            logger.info(f'Skipping {mol.name}.')
            return mol

        try:
            logger.info(f'Optimizing {mol.name}.')
            self.__wrapped__(mol)

        except Exception as ex:
            errormsg = (f'Optimization function '
                        f'"{self.__wrapped__.func.__name__}()" '
                        f'failed on molecule "{mol.name}".')
            logger.error(errormsg, exc_info=True)

        finally:
            mol.optimized = True
            return mol


def do_not_optimize(mol):
    """
    Skips the optimization step.

    This is very useful when debugging so you do not waste your time
    waiting for molecules to get optimized. Use this in the input file
    in place of an optimization function when necessary.

    Parameters
    ----------
    mol : MacroMolecule
        A macromolecule which will not be optimized.

    Returns
    -------
    None : NoneType

    """

    return


def partial_raiser(mol, ofunc):
    """
    Raises and optimizes at random.

    Parameters
    ----------
    mol : MacroMolecule
        The macromolecule being optimized.

    ofunc : FunctionData
        A FunctionData object representing the optimization function
        to be used.

    Returns
    -------
    None : NoneType

    """

    if not np.random.choice([0, 1]):
        raise Exception('Partial raiser.')

    globals()[ofunc.name](mol, **ofunc.params)


def raiser(mol, param1, param2=2):
    """
    Doens't optimize, raises an error instead.

    This function is used to test that when optimization functions
    raise errors during multiprocessing, they are handled correctly.

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

    raise Exception('Raiser optimization function used.')


def rdkit_optimization(mol, embed=False, conformer=-1):
    """
    Optimizes the structure of the molecule using rdkit.

    Uses the ``MMFF`` forcefield.

    Parameters
    ----------
    mol : MacroMolecule
        The macromolecule who's structure should be optimized.

    embed : bool (default = False)
        When ``True`` the structure is guessed before an optimization
        is carried out. This guess structure overrides any previous
        structure.

    conformer : :class:`int`, optional
        The conformer to use.

    Returns
    -------
    None : NoneType

    """

    if conformer == -1:
        conformer = mol.mol.GetConformer(conformer).GetId()

    if embed:
        conf_id = rdkit.EmbedMolecule(mol.mol, clearConfs=False)
        new_conf = rdkit.Conformer(mol.mol.GetConformer(conf_id))
        mol.mol.RemoveConformer(conf_id)
        mol.mol.RemoveConformer(conformer)
        new_conf.SetId(conformer)
        mol.mol.AddConformer(new_conf)

    # Sanitize then optimize the rdkit molecule.
    rdkit.SanitizeMol(mol.mol)
    rdkit.MMFFOptimizeMolecule(mol.mol, confId=conformer)


def rdkit_ETKDG(mol, conformer=-1):
    """
    Does a conformer search with :func:`rdkit.ETKDG` [#]_.

    Parameters
    ----------
    mol : MacroMolecule
        The macromolecule who's structure should be optimized.

    conformer : :class:`int`, optional
        The conformer to use.

    Returns
    -------
    None : NoneType

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    if conformer == -1:
        conformer = mol.mol.GetConformer(conformer).GetId()

    conf_id = rdkit.EmbedMolecule(mol.mol, rdkit.ETKDG())
    new_conf = rdkit.Conformer(mol.mol.GetConformer(conf_id))
    mol.mol.RemoveConformer(conf_id)
    mol.mol.RemoveConformer(conformer)
    new_conf.SetId(conformer)
    mol.mol.AddConformer(new_conf)
