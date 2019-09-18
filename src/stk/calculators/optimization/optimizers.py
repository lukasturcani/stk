"""
Optimizers
==========

#. :class:`.NullOptimizer`
#. :class:`.MMFF`
#. :class:`.UFF`
#. :class:`.ETKDG`
#. :class:`.XTB`
#. :class:`.MacroModelForceField`
#. :class:`.MacroModelMD`
#. :class:`.MOPAC`
#. :class:`.OptimizerSequence`
#. :class:`.CageOptimizerSequence`
#. :class:`.TryCatchOptimizer`
#. :class:`.RaisingOptimizer`


Optimizers are objects used to optimize molecules. Each optimizer is
initialized with some settings and can optimize a molecule
with :meth:`~.Optimizer.optimize`.

.. code-block:: python

    import stk

    mol = stk.BuildingBlock('NCCCN', ['amine'])
    mmff = stk.MMFF()
    mmff.optimize(mol)

    # Optimizers also work with ConstructedMolecule objects.
    polymer = stk.ConstructedMolecule(
        building_blocks=[mol],
        topology_graph=stk.polymer.Linear('A', [0], n=3)
    )
    etkdg = stk.ETKDG()
    etkdg.optimize(polymer)

Sometimes it is desirable to chain multiple optimizations, one after
another. For example, before running an optimization, it may be
desirable to embed a molecule first, to generate an initial structure.
:class:`.OptimizerSequence` may be used for this.

.. code-block:: python

    # Create a new optimizer which chains the previously defined
    # mmff and etkdg optimizers.
    optimizer_sequence = stk.OptimizerSequence(etkdg, mmff)

    # Run each optimizer in sequence.
    optimizer_sequence.optimize(polymer)

By default, running :meth:`.Optimizer.optimize` twice on the same
molecule will perform an optimization a second time on a molecule. If
we want to skip optimizations on molecules which have already been
optimized we can use the `use_cache` flag.

.. code-block:: python

    caching_etkdg = stk.ETKDG(use_cache=True)
    # First optimize call runs an optimization.
    caching_etkdg.optimize(polymer)
    # Second call does nothing.
    caching_etkdg.optimize(polymer)

Caching is done on a per :class:`.Optimizer` basis. Just because the
molecule has been cached by one :class:`.Optimizer` instance does not
mean that a different :class:`.Optimizer` instance will no longer
optimize the molecule.

.. _`adding optimizers`:

Making New Optimizers
---------------------

New optimizers can be made by simply making a class which inherits the
:class:`.Optimizer` class. In addition to this, the new class must
define a :meth:`~.Optimizer.optimize` method. The method must take 1
mandatory `mol` parameter. :meth:`~.Optimizer.optimize` will take the
`mol` and change its structure in whatever way it likes. Beyond this
there are no requirements. New optimizers can be added into the
:mod:`.optimizers` submodule or into a new submodule.

"""

import logging
import numpy as np
import rdkit.Chem.AllChem as rdkit
import warnings
import os
from functools import wraps
import subprocess as sp
import uuid
import shutil
from ...utilities import (
    is_valid_xtb_solvent,
    XTBInvalidSolventError,
    XTBExtractor
)
import pywindow

logger = logging.getLogger(__name__)


def _add_cache_use(optimize):
    """
    Make :meth:`~Optimizer.optimize` use the :attr:`~Optimizer._cache`.

    Decorates `optimize` so that before running it checks if the
    :class:`.Molecule` has already been optimized by the
    optimizer. If so, and :attr:`~Optimizer.use_cache` is ``True``,
    then the molecule is skipped and no optimization is performed.

    Parameters
    ----------
    optimize : :class:`function`
        A function which is to have skipping added to it.

    Returns
    -------
    :class:`function`
        The decorated function.

    """

    @wraps(optimize)
    def inner(self, mol):
        if self._use_cache and mol in self._cache:
            logger.info(f'Skipping optimization on {mol}.')
        else:
            optimize(self, mol)
            if self._use_cache:
                self._cache.add(mol)

    return inner


class Optimizer:
    """
    A base class for optimizers.

    """

    def __init__(self, use_cache=False):
        """
        Initialize an :class:`Optimizer`.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        # Holds every previously optimized molecule if use_cache is
        # true.
        self._cache = set()
        self._use_cache = use_cache

    def __init_subclass__(cls, **kwargs):
        cls.optimize = _add_cache_use(cls.optimize)
        return super().__init_subclass__(**kwargs)

    def is_caching(self):
        """
        ``True`` if the optimizer has caching turned on.

        Returns
        -------
        :class:`bool`
            ``True`` if the optimizer has caching turned on.

        """

        return self._use_cache

    def add_to_cache(self, mol):
        """
        Add a molecule to the cache.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be added to the cache.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._cache.add(mol)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()


class OptimizerSequence(Optimizer):
    """
    Applies optimizers in sequence.

    Examples
    --------
    Let's say we want to embed a molecule with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        optimizer = stk.OptimizerSequence(stk.ETKDG(), stk.MMFF())
        optimizer.optimize(mol)

    """

    def __init__(self, *optimizers, use_cache=False):
        """
        Initialize a :class:`OptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._optimizers = optimizers
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        for optimizer in self._optimizers:
            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol}".')
            optimizer.optimize(mol)


class CageOptimizerSequence(Optimizer):
    """
    Applies :class:`Optimizer` objects to a cage.

    Before each :class:`Optimizer` in the sequence is applied to the
    cage, it is checked to see if it is collapsed. If it is
    collapsed, the optimization sequence ends immediately.

    Examples
    --------
    Let's say we want to embed a cage with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCNCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
        cage = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.cage.FourPlusSix()
        )
        optimizer = stk.CageOptimizerSequence(stk.ETKDG(), stk.MMFF())
        optimizer.optimize(cage)

    """

    def __init__(self, *optimizers, use_cache=False):
        """
        Initialize a :class:`CageOptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            The :class:`Optimizers` used in sequence to optimize
            cage molecules.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._optimizers = optimizers
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The cage to be optimized. It must have an attribute
            :attr:`num_windows`, which contains the expected number
            of windows if the cage is not collapsed.

        Returns
        -------
        None : :class:`NoneType`

        """

        for optimizer in self._optimizers:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                loader = pywindow.molecular.Molecule.load_rdkit_mol
                pw_molecule = loader(mol.to_rdkit_mol())
                windows = pw_molecule.calculate_windows()
            logger.debug(f'Windows found: {windows}.')

            if windows is None or len(windows) != mol.num_windows:
                logger.info(f'"{mol}" is collapsed, exiting early.')
                return

            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol}".')
            optimizer.optimize(mol)


class NullOptimizer(Optimizer):
    """
    Does not perform optimizations.

    """

    def optimize(self, mol):
        """
        Do not optimize `mol`.

        This function just returns immediately without changing the
        molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule.

        Returns
        -------
        None : :class:`NoneType`

        """

        return


class TryCatchOptimizer(Optimizer):
    """
    Try to optimize with a :class:`Optimizer`, use another on failure.

    Examples
    --------
    .. code-block:: python

        import stk

        # Create some molecules to optimize.
        mol1 = stk.BuildingBlock('NCCN', ['amine'])
        mol2 = stk.BuildingBlock('CCCCC')
        mol3 = stk.BuildingBlock('O=CCCN')

        # Create an optimizer which may fail.
        uff = stk.UFF()

        # Create a backup optimizer.
        mmff = stk.MMFF()

        # Make an optimizer which tries to run raiser and if that
        # raises an error, will run mmff on the molecule instead.
        try_catch = stk.TryCatchOptimizer(
            try_optimizer=uff,
            catch_optimizer=mmff
        )

        # Optimize the molecules. In each case if the optimization with
        # UFF fails, MMFF is used to optimize the molecule instead.
        try_catch.optimize(mol1)
        try_catch.optimize(mol2)
        try_catch.optimzie(mol3)

    """

    def __init__(
        self,
        try_optimizer,
        catch_optimizer,
        use_cache=False
    ):
        """
        Initialize a :class:`TryCatchOptimizer` instance.

        Parameters
        ----------
        try_optimizer : :class:`Optimizer`
            The optimizer which is used initially to try and optimize a
            :class:`.Molecule`.

        catch_optimizer : :class:`Optimizer`
            If :attr:`try_optimizer` raises an error, this optimizer is
            run on the :class:`.Molecule` instead.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._try_optimizer = try_optimizer
        self._catch_optimizer = catch_optimizer
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        try:
            return self._try_optimizer.optimize(mol)
        except Exception:
            try_name = self._try_optimizer.__class__.__name__
            catch_name = self._catch_optimizer.__class__.__name__
            logger.error(
                f'{try_name} failed, trying {catch_name}.',
                exc_info=True
            )
            return self._catch_optimizer.optimize(mol)


class RaisingOptimizerError(Exception):
    ...


class RaisingOptimizer(Optimizer):
    """
    Raises and optimizes at random.

    This optimizer is used for debugging to simulate optimizations
    which sometimes complete successfully and sometimes
    randomly fail.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        partial_raiser = stk.RaisingOptimizer(
            optimizer=stk.ETKDG(),
            fail_chance=0.75
        )
        # 75 % chance an error will be raised by calling optimize.
        partial_raiser.optimize(mol)

    """

    def __init__(self, optimizer, fail_chance=0.5, use_cache=False):
        """
        Initialize :class:`PartialRaiser`.

        Parameters
        ----------
        optimizer : :class:`Optimizer`
            When the optimizer does not fail, it uses this
            :class:`Optimizer` to optimize molecules.

        fail_chance : :class:`float`, optional
            The probability that the optimizer will raise an error each
            time :meth:`optimize` is used.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._optimizer = optimizer
        self._fail_chance = fail_chance
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`RaisingOptimizerError`
            This error is raised at random.

        """

        if np.random.rand() < self._fail_chance:
            raise RaisingOptimizerError('Used RaisingOptimizer.')
        return self._optimizer.optimize(mol)


class MMFF(Optimizer):
    """
    Use the MMFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        mmff = stk.MMFF()
        mmff.optimize(mol)

    """

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.MMFFOptimizeMolecule(rdkit_mol)
        mol.update_from_rdkit_mol(rdkit_mol)


class UFF(Optimizer):
    """
    Use the UFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        uff = stk.UFF()
        uff.optimize(mol)

    """

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        rdkit_mol = mol.to_rdkit_mol()
        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(rdkit_mol)
        rdkit.UFFOptimizeMolecule(rdkit_mol)
        mol.update_from_rdkit_mol(rdkit_mol)


class ETKDG(Optimizer):
    """
    Uses the ETKDG [#]_ v2 algorithm to find an optimized structure.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        etkdg = stk.ETKDG()
        etkdg.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, random_seed=12, use_cache=False):
        """
        Initialize a :class:`ETKDG` instance.

        Parameters
        ----------
        random_seed : :class:`int`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        self._random_seed = random_seed
        super().__init__(use_cache=use_cache)

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        params = rdkit.ETKDGv2()
        params.clearConfs = True
        params.random_seed = self._random_seed

        rdkit_mol = mol.to_rdkit_mol()
        rdkit.EmbedMolecule(rdkit_mol, params)
        mol.set_position_matrix(
            position_matrix=rdkit_mol.GetConformer().GetPositions()
        )


class XTBOptimizerError(Exception):
    ...


class XTBConvergenceError(XTBOptimizerError):
    ...


class XTB(Optimizer):
    """
    Uses GFN-xTB [1]_ to optimize molecules.

    Notes
    -----
    When running :meth:`optimize`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised, so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as :class:`.XTB` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Furthermore, :meth:`optimize` will check that the
    structure is adequately optimized by checking for negative
    frequencies after a Hessian calculation. `max_runs` can be
    provided to the initializer to set the maximum number of
    optimizations which will be attempted at the given
    `opt_level` to obtain an optimized structure. However, we outline
    in the examples how to iterate over `opt_levels` to increase
    convergence criteria and hopefully obtain an optimized structure.
    The presence of negative frequencies can occur even when the
    optimization has converged based on the given `opt_level`.

    Attributes
    ----------
    incomplete : :class:`set` of :class:`.Molecule`
        A :class:`set` of molecules passed to :meth:`optimize` whose
        optimzation was incomplete.

    Examples
    --------
    Note that for :class:`.ConstructedMolecule` objects constructed by
    ``stk``, :class:`XTB` should usually be used in a
    :class:`.OptimizerSequence`. This is because xTB only uses
    xyz coordinates as input and so will not recognize the long bonds
    created during construction. An optimizer which can minimize
    these bonds should be used before :class:`XTB`.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCNCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CCCC=O', ['aldehyde'])
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear("AB", [0, 0], 3)
        )

        xtb = stk.OptimizerSequence(
            stk.UFF(),
            stk.XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
        )
        xtb.optimize(polymer)

    By default, all optimizations with xTB are performed using the
    ``--ohess`` flag, which forces the calculation of a numerical
    Hessian, thermodynamic properties and vibrational frequencies.
    :meth:`optimize` will check that the structure is appropriately
    optimized (i.e. convergence is obtained and no negative vibrational
    frequencies are present) and continue optimizing a structure (up to
    `max_runs` times) until this is achieved. This loop, by
    default, will be performed at the same `opt_level`. The
    following example shows how a user may optimize structures with
    tigher convergence criteria (i.e. different `opt_level`)
    until the structure is sufficiently optimized. Furthermore, the
    calculation of the Hessian can be turned off using
    `max_runs` to ``1`` and `calculate_hessian` to ``False``.

    .. code-block:: python

        # Use crude optimization with max_runs=1 because this will
        # not achieve optimization and rerunning it is unproductive.
        xtb_crude = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_crude',
            unlimited_memory=True,
            opt_level='crude',
            max_runs=1,
            calculate_hessian=True
        )
        # Use normal optimization with max_runs == 2.
        xtb_normal = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_normal',
            unlimited_memory=True,
            opt_level='normal',
            max_runs=2
        )
        # Use vtight optimization with max_runs == 2, which should
        # achieve sufficient optimization.
        xtb_vtight = stk.XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_vtight',
            unlimited_memory=True,
            opt_level='vtight',
            max_runs=2
        )

        optimizers = [xtb_crude, xtb_normal, xtb_vtight]
        for optimizer in optimizers:
            optimizer.optimize(polymer)
            if polymer not in optimizer.incomplete:
                break

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/setup.html

    """

    def __init__(
        self,
        xtb_path,
        gfn_version=2,
        output_dir=None,
        opt_level='normal',
        max_runs=2,
        calculate_hessian=True,
        num_cores=1,
        electronic_temperature=300,
        solvent=None,
        solvent_grid='normal',
        charge=0,
        num_unpaired_electrons=0,
        unlimited_memory=False,
        use_cache=False
    ):
        """
        Initialize a :class:`XTB` instance.

        Parameters
        ----------
        xtb_path : :class:`str`
            The path to the xTB executable.

        gfn_version : :class:`int`, optional
            Parameterization of GFN to use in xTB.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/basics.html.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
            ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/optimization.html
            .

        max_runs : :class:`int`, optional
            Maximum number of optimizations to attempt in a row.

        calculate_hessian : :class:`bool`, optional
            Toggle calculation of the hessian and vibrational
            frequencies after optimization. ``True`` is required to
            check that the structure is completely optimized.
            ``False`` will drastically speed up the calculation but
            potentially provide incomplete optimizations and forces
            :attr:`max_runs` to be ``1``.

        num_cores : :class:`int`, optional
            The number of cores xTB should use.

        electronic_temperature : :class:`int`, optional
            Electronic temperature in Kelvin.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit
            solvent.
            Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
            or ``'extreme'``.
            For details see
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

        charge : :class:`int`, optional
            Formal molecular charge.

        num_unpaired_electrons : :class:`int`, optional
            Number of unpaired electrons.

        unlimited_memory : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without
            constraints on the stack size. If memory issues are
            encountered, this should be ``True``, however this may
            raise issues on clusters.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule.

        """

        if solvent is not None:
            solvent = solvent.lower()
            if gfn_version == 0:
                raise XTBInvalidSolventError(
                    f'No solvent valid for version',
                    f' {gfn_version!r}.'
                )
            if not is_valid_xtb_solvent(gfn_version, solvent):
                raise XTBInvalidSolventError(
                    f'Solvent {solvent!r} is invalid for ',
                    f'version {gfn_version!r}.'
                )

        if not calculate_hessian and max_runs != 1:
            max_runs = 1
            logger.warning(
                'Requested that hessian calculation was skipped '
                'but the number of optimizations requested was '
                'greater than 1. The number of optimizations has been '
                'set to 1.'
            )

        self._xtb_path = xtb_path
        self._gfn_version = str(gfn_version)
        self._output_dir = output_dir
        self._opt_level = opt_level
        self._max_runs = max_runs
        self._calculate_hessian = calculate_hessian
        self._num_cores = str(num_cores)
        self._electronic_temperature = str(electronic_temperature)
        self._solvent = solvent
        self._solvent_grid = solvent_grid
        self._charge = str(charge)
        self._num_unpaired_electrons = str(num_unpaired_electrons)
        self._unlimited_memory = unlimited_memory
        self.incomplete = set()
        super().__init__(use_cache=use_cache)

    def _has_neg_frequencies(self, output_file):
        """
        Check for negative frequencies.

        Parameters
        ----------
        output_file : :class:`str`
            Name of output file with xTB results.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if a negative frequency is present.

        """
        xtbext = XTBExtractor(output_file=output_file)
        # Check for one negative frequency, excluding the first
        # 6 frequencies.
        return any(x < 0 for x in xtbext.frequencies[6:])

    def _is_complete(self, output_file):
        """
        Check if xTB optimization has completed and converged.

        Parameters
        ----------
        output_file : :class:`str`
            Name of xTB output file.

        Returns
        -------
        :class:`bool`
            Returns ``False`` if a negative frequency is present.

        Raises
        -------
        :class:`XTBOptimizerError`
            If the optimization failed.

        :class:`XTBConvergenceError`
            If the optimization did not converge.

        """
        if output_file is None:
            # No simulation has been run.
            return False
        # If convergence is achieved, then .xtboptok should exist.
        if os.path.exists('.xtboptok'):
            # Check for negative frequencies in output file if the
            # hessian was calculated.
            # Return True if there exists at least one.
            if self._calculate_hessian:
                return not self._has_neg_frequencies(output_file)
            else:
                return True
        elif os.path.exists('NOT_CONVERGED'):
            raise XTBConvergenceError('Optimization not converged.')
        else:
            raise XTBOptimizerError('Optimization failed to complete')

    def _run_xtb(self, xyz, out_file):
        """
        Run GFN-xTB.

        Parameters
        ----------
        xyz : :class:`str`
            The name of the input structure ``.xyz`` file.

        out_file : :class:`str`
            The name of output file with xTB results.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Modify the memory limit.
        if self._unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        # Set optimization level and type.
        if self._calculate_hessian:
            # Do optimization and check hessian.
            optimization = f'--ohess {self._opt_level}'
        else:
            # Do optimization.
            optimization = f'--opt {self._opt_level}'

        if self._solvent is not None:
            solvent = f'--gbsa {self._solvent} {self._solvent_grid}'
        else:
            solvent = ''

        cmd = (
            f'{memory} {self._xtb_path} {xyz} '
            f'--gfn {self._gfn_version} '
            f'{optimization} --parallel {self._num_cores} '
            f'--etemp {self._electronic_temperature} '
            f'{solvent} --chrg {self._charge} '
            f'--uhf {self._num_unpaired_electrons}'
        )

        with open(out_file, 'w') as f:
            # Note that sp.call will hold the program until completion
            # of the calculation.
            sp.call(
                cmd,
                stdin=sp.PIPE,
                stdout=f,
                stderr=sp.PIPE,
                # Shell is required to run complex arguments.
                shell=True
            )

    def _run_optimizations(self, mol):
        """
        Run loop of optimizations on `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """
        for run in range(self._max_runs):
            xyz = f'input_structure_{run+1}.xyz'
            out_file = f'optimization_{run+1}.output'
            mol.write(xyz)
            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            coord_file = 'xtbhess.coord'
            coord_exists = os.path.exists(coord_file)
            output_xyz = 'xtbopt.xyz'
            opt_complete = self._is_complete(out_file)
            if not opt_complete:
                if coord_exists:
                    # The calculation is incomplete.
                    # Update mol from xtbhess.coord and continue.
                    mol.update_from_file(coord_file)
                else:
                    # Update mol from xtbopt.xyz.
                    mol.update_from_file(output_xyz)
                    # If the negative frequencies are small, then GFN
                    # may not produce the restart file. If that is the
                    # case, exit optimization loop and warn.
                    self.incomplete.add(mol)
                    logging.warning(
                        f'Small negative frequencies present in {mol}.'
                    )
                    return False
            else:
                # Optimization is complete.
                # Update mol from xtbopt.xyz.
                mol.update_from_file(output_xyz)
                break
        return opt_complete

    def optimize(self, mol):
        """
        Optimize `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Remove mol from self.incomplete if present.
        if mol in self.incomplete:
            self.incomplete.remove(mol)

        if self._output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self._output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            complete = self._run_optimizations(mol)
        finally:
            os.chdir(init_dir)

        if not complete:
            self.incomplete.add(mol)
            logging.warning(f'Optimization is incomplete for {mol}.')
