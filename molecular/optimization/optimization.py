"""
Defines optimization functions.

Extending MMEA: Adding optimization functions
---------------------------------------------
New optimization functions are added by writing them into this module.
The only requirement is that the first argument is ``macro_mol`` and
there is a keyword argument ``logger``. This allows users to identify
which arguments are handled automatically by MMEA and which need to be
defined in the input file. The convention is that if the optimization
function takes the arguments ``macro_mol`` and ``logger``, the user
does not have to specify them argument in the input file.

When defining optimization functions the ``logger`` argument should be
used for logging as a normal logger from the ``logging`` library would.
When running the GA a special logger compatible with multiprocessing
is automatically placed in this argument. It may be useful to define
the logger argument as a keyword argument::

    opt_func(macro_mol, somearg, logger=logging.getLogger(__name__)):
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
from traceback import format_exc
import psutil

from .macromodel import macromodel_opt, macromodel_cage_opt

from ...convenience_tools import FakeLogger, mplogger, StopLogging


logger = logging.getLogger(__name__)


def _optimize_all(func_data, population):
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

    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the
    # module.
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may
    # require.
    p_func = _OptimizationFunc(partial(func, **func_data.params),
                               logq=logq)

    # Apply the function to every member of the population, in
    # parallel.
    with mp.get_context('spawn').Pool(psutil.cpu_count(False)) as pool:
        optimized = pool.map(p_func, population)
    # Make sure the cache is updated with the optimized versions.
    for member in optimized:
        member.update_cache()

    logq.put((exit_, exit_))
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

    def __init__(self, func, logq=None):
        wraps(func)(self)
        self.logq = logq

    def __call__(self, macro_mol, *args, **kwargs):
        logger = logging.getLogger(__name__)
        if self.logq is not None:
            logger = FakeLogger(self.logq)

        if macro_mol.optimized:
            logger.info('Skipping {}.'.format(macro_mol.name))
            return macro_mol

        try:
            logger.info('Optimizing {}.'.format(macro_mol.name))
            self.__wrapped__(macro_mol, *args, **kwargs, logger=logger)

        except Exception as ex:
            errormsg = ('Optimization function "{}()" '
                        'failed on molecule "{}".').format(
                        self.__wrapped__.func.__name__, macro_mol.name)
            logger.error((errormsg+'\n'+format_exc()).strip())

        finally:
            macro_mol.optimized = True
            return macro_mol


def do_not_optimize(macro_mol, logger=logger):
    """
    Skips the optimization step.

    This is very useful when debugging so you do not waste your time
    waiting for molecules to get optimized. Use this in the input file
    in place of an optimization function when necessary.

    Parameters
    ----------
    macro_mol : MacroMolecule
        A macromolecule which will not be optimized.

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Returns
    -------
    None : NoneType

    """

    return


def partial_raiser(macro_mol, ofunc, logger=logger):
    """
    Raises and optimizes at random.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule being optimized.

    ofunc : FunctionData
        A FunctionData object representing the optimization function
        to be used.

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Modifies
    --------
    macro_mol.mol
        The rdkit instance is optimized.

    Returns
    -------
    None : NoneType

    """

    if not np.random.choice([0, 1]):
        raise Exception('Partial raiser.')

    globals()[ofunc.name](macro_mol, **ofunc.params)


def raiser(macro_mol, param1, param2=2, logger=logger):
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

    raise Exception('Raiser optimization function used.')


def rdkit_optimization(macro_mol, embed=False, logger=logger):
    """
    Optimizes the structure of the molecule using rdkit.

    Uses the MMFF forcefield.

    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure should be optimized.

    embed : bool (default = False)
        When ``True`` the structure is guessed before an optimization
        is carried out. This guess structure overrides any previous
        structure.

    logger : FakeLogger or logging.Logger, optional
        Used for logging. Not used by this function.

    Modifies
    --------
    macro_mol.mol
        The rdkit molecule held in this attribute has it's structure
        changed as a result of the optimization. This means the
        ``Conformer`` instance held by the rdkit molecule is changed.

    Returns
    -------
    None : NoneType

    """

    if embed:
        rdkit.EmbedMolecule(macro_mol.mol)

    # Sanitize then optimize the rdkit molecule.
    rdkit.SanitizeMol(macro_mol.mol)
    rdkit.MMFFOptimizeMolecule(macro_mol.mol)
