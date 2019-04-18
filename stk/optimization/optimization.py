"""
Defines optimization functions.

.. _`adding optimization functions`:

Extending stk: Adding optimization functions.
---------------------------------------------

New optimization functions are added by writing them in this module.
The only requirement is that the first argument is `mol`. This allows
users to identify which arguments are handled automatically by ``stk``
and which need to be defined in the input file. The convention is that
if the optimization function takes the argument `mol`, the user does
not have to specify it in the input file.

An optimization function should update the ``rdkit`` molecule in
:attr:`.Molecule.mol`. The return values of optimization functions
are discarded by the GA.

Optimizations can be complicated. If the use of helper functions is
required make sure that they are private, i.e. that their names begin
with a leading underscore. In the event that the optimization is so
complex that it requires its own module or file, place it in the same
folder as this file. Then import the optimization function here. Make
sure that only the optimization functions are imported back into this
file, not any of the helper functions or classes.

"""

import rdkit.Chem.AllChem as rdkit
import multiprocessing as mp
from functools import partial, wraps
import numpy as np
import logging
from threading import Thread

from .macromodel import macromodel_opt, macromodel_cage_opt
from ..utilities import daemon_logger, logged_call, OPTIONS


logger = logging.getLogger(__name__)


def _optimize_all(fn, population, processes):
    """
    Run opt function on all population members in parallel.

    Parameters
    ----------
    fn : :class:`function`
        An optimization function, which takes a single argument,
        the molecule to be optimized.

    population : :class:`.Population`
        The :class:`.Population` instance who's members are to be
        optimized.

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

    opt_fn = _OptimizationFunc(fn)

    # Apply the function to every member of the population, in
    # parallel.
    with mp.get_context('spawn').Pool(processes) as pool:
        optimized = pool.starmap(logged_call,
                                 ((logq, opt_fn, mem) for
                                  mem in population))

    # Update the structures in the population.
    sorted_opt = sorted(optimized, key=lambda m: m.key)
    sorted_pop = sorted(population, key=lambda m: m.key)
    for old, new in zip(sorted_pop, sorted_opt):
        old.__dict__ = dict(vars(new))

    # Make sure the cache is updated with the optimized versions.
    if OPTIONS['cache']:
        for member in optimized:
            member.update_cache()

    logq.put(None)
    log_thread.join()


def _optimize_all_serial(fn, population):
    """
    Run opt function on all population members sequentially.

    Parameters
    ----------
    fn : :class:`function`
        An optimization function, which takes a single argument,
        the molecule to be optimized.

    population : :class:`.Population`
        The :class:`.Population` instance who's members are to be
        optimized.

    Returns
    -------
    None : :class:`NoneType`

    """

    opt_fn = _OptimizationFunc(fn)

    # Apply the function to every member of the population.
    for member in population:
        opt_fn(member)


class _OptimizationFunc:
    """
    A decorator for optimization functions.

    This decorator is applied to all optimization functions
    automatically in :func:`_optimize_all`. It should not be applied
    explicitly when defining the functions.

    This decorator prevents optimization functions from raising if
    they fail (necessary for multiprocessing) and prevents them from
    being run twice on the same molecule.

    """

    def __init__(self, func):
        wraps(func)(self)

    def __call__(self, mol):
        """
        Decorates and calls the optimization function.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        :class:`.Molecule`
            The optimized molecule.

        """

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


def partial_raiser(mol, fn):
    """
    Raises and optimizes at random.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule being optimized.

    fn : :class:`function`
        An optimization function to be used, must take a single
        argument, the molecule to be optimized.

    Returns
    -------
    None : :class:`.NoneType`

    """

    if not np.random.choice([0, 1]):
        raise Exception('Partial raiser.')

    fn(mol)


def raiser(mol, param1, param2=2):
    """
    Doens't optimize, raises an error instead.

    This function is used to test that when optimization functions
    raise errors during multiprocessing, they are handled correctly.

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

    raise Exception('Raiser optimization function used.')


def rdkit_optimization(mol, embed=False, conformer=-1):
    """
    Optimizes the structure of the molecule using ``rdkit``.

    Uses the ``MMFF`` forcefield.

    Parameters
    ----------
    mol : :class:`.Molecule`
        The molecule who's structure should be optimized.

    embed : :class:`bool`, optional
        When ``True`` the structure is guessed before an optimization
        is carried out. This guess structure overrides any previous
        structure.

    conformer : :class:`int`, optional
        The conformer to use.

    Returns
    -------
    None : :class:`NoneType`

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
    mol : :class:`.Molecule`
        The molecule who's structure should be optimized.

    conformer : :class:`int`, optional
        The conformer to use.

    Returns
    -------
    None : :class:`NoneType`

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
