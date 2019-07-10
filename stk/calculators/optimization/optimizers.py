"""
Defines optimizers.

See :mod:`.optimization`.

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
from ...utilities import (is_valid_xtb_solvent,
                          XTBInvalidSolventError,
                          XTBExtractor)

logger = logging.getLogger(__name__)


def _add_cache_use(optimize):
    """
    Makes :meth:`~Optimizer.optimize` use the :attr:`~Optimizer.cache`.

    Decorates `optimize` so that before running it checks if the
    :class:`.Molecule` and conformer have already been optimized by the
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
    def inner(self, mol, conformer=-1):
        key = (mol.key, conformer)
        if self.use_cache and key in self.cache:
            logger.info(
                f'Skipping optimization on '
                f'"{mol.name}" conformer {conformer}.'
            )
        else:
            optimize(self, mol, conformer)
            if self.use_cache:
                self.cache.add(key)

    return inner


class Optimizer:
    """
    A base class for optimizers.

    Attributes
    ----------
    cache : :class:`set`
        A :class:`set` of the form

        .. code-block:: python

            cache = {
                (mol1, conformer1),
                (mol1, conformer2),
                (mol2, conformer2)
            }

        which holds every :class:`Molecule` and conformer optimized
        by the :class:`Optimizer`. Here ``mol1`` and ``mol2`` are
        :class:`.Molecule` objects and ``conformer1`` and
        ``conformer2`` are :class:`int`, which are the ids of the
        optimized conformers of the molecules.

    use_cache : :class:`bool`
        If ``True`` :meth:`optimize` will not run twice on the same
        molecule and conformer.

    """

    def __init__(self, use_cache=False):
        """
        Initializes an :class:`Optimizer`.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        """

        self.cache = set()
        self.use_cache = use_cache

    def __init_subclass__(cls, **kwargs):
        cls.optimize = _add_cache_use(cls.optimize)
        return super().__init_subclass__(**kwargs)

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()


class OptimizerSequence(Optimizer):
    """
    Applies optimizers in sequence.

    Attributes
    ----------
    optimizers : :class:`tuple` of :class:`Optimizer`
        A number of optimizers, each of which gets applied to a
        molecule, based on the order in this :class:`tuple`.

    Examples
    --------
    Let's say we want to embed a molecule with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit
        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        etkdg = RDKitEmbedder(rdkit.ETKDG())
        mmff = RDKitForceField(rdkit.MMFFOptimizeMolecule)
        optimizer = OptimizerSequence(etkdg, mmff)
        optimizer.optimize(mol)

    """

    def __init__(self, *optimizers, use_cache=False):
        """
        Initializes a :class:`OptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        """

        self.optimizers = optimizers
        super().__init__(use_cache=use_cache)

    def optimize(self, mol, conformer=-1):
        """
        Chains multiple :class:`Optimizer` instances together.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        for optimizer in self.optimizers:
            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol.name}".')
            optimizer.optimize(mol, conformer)


class CageOptimizerSequence(Optimizer):
    """
    Applies :class:`Optimizer` objects to a :class:`.Cage`.

    Before each :class:`Optimizer` in the sequence is applied to the
    :class:`.Cage`, it is checked to see if it is collapsed. If it is
    collapsed, the optimization sequence ends immediately.

    Attributes
    ----------
    optimizers : :class:`tuple` of :class:`Optimizer`
        The :class:`Optimizer` objects which are used to optimize a
        :class:`.Cage` molecule.

    Examples
    --------
    Let's say we want to embed a cage with ETKDG first and then
    minimize it with the MMFF force field.

    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit
        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        etkdg = RDKitEmbedder(rdkit.ETKDG())
        mmff = RDKitForceField(rdkit.MMFFOptimizeMolecule)
        optimizer = CageOptimizerSequence(etkdg, mmff)
        optimizer.optimize(mol)

    """

    def __init__(self, *optimizers, use_cache=False):
        """
        Initializes a :class:`CageOptimizerSequence` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            The :class:`Optimizers` used in sequence to optimize
            :class:`.Cage` molecules.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        """

        self.optimizers = optimizers
        super().__init__(use_cache=use_cache)

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Cage`
            The cage to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        for optimizer in self.optimizers:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                windows = mol.windows(conformer)

            logger.debug(f'Windows found: {windows}.')
            expected_windows = mol.topology.n_windows
            if windows is None or len(windows) != expected_windows:
                logger.info(
                    f'"{mol.name}" is collapsed, exiting early.'
                )
                return

            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol.name}".')
            optimizer.optimize(mol, conformer)


class NullOptimizer(Optimizer):
    """
    Does not perform optimizations.

    """

    def optimize(self, mol, conformer=-1):
        """
        Does not optimize a molecule.

        This function just returns immediately without changing the
        molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            A molecule.

        conformer : :class:`int`, optional
            A conformer.

        Returns
        -------
        None : :class:`NoneType`

        """

        return


class TryCatchOptimizer(Optimizer):
    """
    Try to optimize with a :class:`Optimizer`, use another on failure.

    Attributes
    ----------
    try_optimizer : :class:`Optimizer`
        The optimizer which is used initially to try and optimize a
        :class:`.Molecule`.

    catch_optimizer : :class:`Optimizer`
        If :attr:`try_optimizer` raises an error, this optimizer is
        run on the :class:`.Molecule` instead.

    Examples
    --------
    .. code-block:: python

        # Create some molecules to optimize.
        mol1 = StructUnit2.smiles_init('NCCN', []'amine'])
        mol2 = StructUnit.smiles_init('CCCCC')
        mol3 = StructUnit.smiles_init('O=CCCN')

        # Create an optimizer which may fail.
        uff = UFF()

        # Create a backup optimizer.
        mmff = MMFF()

        # Make an optimizer which tries to run raiser and if that
        # raises an error, will run mmff on the molecule instead.
        try_catch = TryCatchOptimizer(try_optimizer=uff,
                                      catch_optimizer=mmff)

        # Optimize the molecules. In each case if the optimization with
        # UFF fails, MMFF is used to optimize the molecule instead.
        try_catch.optimize(mol1)
        try_catch.optimize(mol2)
        try_catch.optimzie(mol3)

    """

    def __init__(self,
                 try_optimizer,
                 catch_optimizer,
                 use_cache=False):
        """
        Initializes a :class:`TryCatchOptimizer` instance.

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
            molecule and conformer.

        """

        self.try_optimizer = try_optimizer
        self.catch_optimizer = catch_optimizer
        super().__init__(use_cache=use_cache)

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        try:
            return self.try_optimizer.optimize(mol, conformer)
        except Exception:
            try_name = self.try_optimizer.__class__.__name__
            catch_name = self.catch_optimizer.__class__.__name__
            logger.error(
                f'{try_name} failed, trying {catch_name}.',
                exc_info=True
            )
            return self.catch_optimizer.optimize(mol, conformer)


class RaisingOptimizerError(Exception):
    ...


class RaisingOptimizer(Optimizer):
    """
    Raises and optimizes at random.

    This optimizer is used for debugging to simulate optimization
    functions which sometimes completes successfully and sometimes
    randomly fails.

    Attributes
    ----------
    optimizer : :class:`Optimizer`
        When the optimizer does not fail, it uses this
        :class:`Optimizer` to optimize molecules.

    fail_chance : :class:`float`
        The probability that the optimizer will raise an error each
        time :meth:`optimize` is used.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        etkdg = RDKitEmbedder(rdkit.ETKDG())
        partial_raiser = RaisingOptimizer(etkdg, fail_chance=0.75)
        # 75 % chance an error will be raised by calling optimize.
        partial_raiser.optimize(mol)

    """

    def __init__(self, optimizer, fail_chance=0.5, use_cache=False):
        """
        Initializes :class:`PartialRaiser`.

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
            molecule and conformer.

        """

        self.optimizer = optimizer
        self.fail_chance = fail_chance
        super().__init__(use_cache=use_cache)

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`RaisingOptimizerError`
            This error is raised at random.

        """

        if np.random.rand() < self.fail_chance:
            raise RaisingOptimizerError('Used RaisingOptimizer.')
        return self.optimizer.optimize(mol)


class MMFF(Optimizer):
    """
    Use the MMFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        mmff = MMFF()
        mmff.optimize(mol)

    """

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule with the MMFF force field.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(mol.mol)
        rdkit.MMFFOptimizeMolecule(mol.mol, confId=conformer)


class UFF(Optimizer):
    """
    Use the UFF force field to optimize molecules.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        uff = UFF()
        uff.optimize(mol)

    """

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule with the UFF force field.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        # Needs to be sanitized to get force field params.
        rdkit.SanitizeMol(mol.mol)
        rdkit.UFFOptimizeMolecule(mol.mol, confId=conformer)


class ETKDG(Optimizer):
    """
    Uses the ETKDG [#]_ v2 algorithm to find an optimized structure.

    Attributes
    ----------
    random_seed : :class:`int`
        The random seed to use.

    Examples
    --------
    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit
        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        etkdg = ETKDG()
        etkdg.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, random_seed=12, use_cache=False):
        """
        Initializes a :class:`ETKDG` instance.

        Parameters
        ----------
        random_seed : :class:`int`, optional
            The random seed to use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        """

        self.random_seed = random_seed
        super().__init__(use_cache=use_cache)

    def optimize(self, mol, conformer=-1):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()

        params = rdkit.ETKDGv2()
        params.clearConfs = False
        params.random_seed = self.random_seed

        conf_id = rdkit.EmbedMolecule(mol.mol, params)
        # Make sure that the conformer order is not re-arranged.
        positions = mol.position_matrix(conformer=conf_id)
        mol.set_position_from_matrix(positions, conformer=conformer)
        mol.mol.RemoveConformer(conf_id)


class XTBOptimizerError(Exception):
    ...


class XTBConvergenceError(XTBOptimizerError):
    ...


class XTB(Optimizer):
    """
    Uses GFN-xTB to optimize molecules.

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
    frequencies after a Hessian calculation. :attr:`max_runs`
    optimizations will be attempted at the given :attr:`opt_level`
    to obtain an optimized structure. However, we outline in the
    examples how to iterate over :attr:`opt_levels` to increase
    convergence criteria and hopefully obtain an optimized structure.
    The presence of negative frequencies can occur even when the
    optimization has converged based on the given
    :attr:`opt_level`.

    Attributes
    ----------
    xtb_path : :class:`str`
        The path to the xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN to use in xTB.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/basics.html.

    output_dir : :class:`str`
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    opt_level : :class:`str`
        Optimization level to use.
        Can be one of ``'crude'``, ``'sloppy'``, ``'loose'``,
        ``'lax'``, ``'normal'``, ``'tight'``, ``'vtight'``
        or ``'extreme'``.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/optimization.html.

    max_runs : :class:`int`
        Maximum number of optimizations to attempt in a row.

    calculate_hessian : :class:`bool`
        Toggle calculation of the hessian and vibrational frequencies
        after optimization. ``True`` is required to check that the
        structure is completely optimized. ``False`` will drastically
        speed up the calculation but potentially provide incomplete
        optimizations and forces :attr:`max_runs` to be ``1``.

    num_cores : :class:`str`
        The number of cores xTB should use.

    electronic_temperature : :class:`str`
        Electronic temperature in Kelvin.

    solvent : :class:`str`
        Solvent to use in GBSA implicit solvation method.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

    solvent_grid : :class:`str`
        Grid level to use in SASA calculations for GBSA implicit
        solvent.
        Can be one of ``'normal'``, ``'tight'``, ``'verytight'``
        or ``'extreme'``.
        For details see
        https://xtb-docs.readthedocs.io/en/latest/gbsa.html.

    charge : :class:`str`
        Formal molecular charge.

    num_unpaired_electrons : :class:`str`
        Number of unpaired electrons.

    unlimited_memory : :class:`bool`
        If ``True`` :meth:`optimize` will be run without constraints on
        the stack size. If memory issues are encountered, this should
        be ``True``, however this may raise issues on clusters.

    incomplete : :class:`set`
        A :class:`set` holding tuples of the form ``(mol, conformer)``,
        which are the :class:`.Molecule` objects and conformer ids
        passed to :meth:`optimize` whose optimzation was incomplete.

    Examples
    --------
    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`XTB` should usually be used in a
    :class:`OptimizerSequence`. This is because xTB only uses
    xyz coordinates as input and so will not recognize the long bonds
    created during assembly. An optimizer which can minimize
    these bonds should be used before :class:`XTB`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        xtb = OptimizerSequence(
            UFF(),
            XTB(xtb_path='/opt/gfnxtb/xtb', unlimited_memory=True)
        )
        xtb.optimize(polymer)

    By default, all optimizations with xTB are performed using the
    ``--ohess`` flag, which forces the calculation of a numerical
    Hessian, thermodynamic properties and vibrational frequencies.
    :meth:`optimize` will check that the structure is appropriately
    optimized (i.e. convergence is obtained and no negative vibrational
    frequencies are present) and continue optimizing a structure (up to
    :attr:`max_runs` times) until this is achieved. This loop, by
    default, will be performed at the same :attr:`opt_level`. The
    following example shows how a user may optimize structures with
    tigher convergence criteria (i.e. different :attr:`opt_level`)
    until the structure is sufficiently optimized. Furthermore, the
    calculation of the Hessian can be turned off using
    :attr:`max_runs` to ``1`` and :attr:`calculate_hessian` to
    ``False``.

    .. code-block:: python

        # Use crude optimization with max_runs=1 because this will
        # not achieve optimization and rerunning it is unproductive.
        xtb_crude = XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_crude',
            unlimited_memory=True,
            opt_level='crude',
            max_runs=1,
            calculate_hessian=True
        )
        # Use normal optimization with max_runs == 2.
        xtb_normal = XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_normal',
            unlimited_memory=True,
            opt_level='normal',
            max_runs=2
        )
        # Use vtight optimization with max_runs == 2, which should
        # achieve sufficient optimization.
        xtb_vtight = XTB(
            xtb_path='/opt/gfnxtb/xtb',
            output_dir='xtb_vtight',
            unlimited_memory=True,
            opt_level='vtight',
            max_runs=2
        )

        # The conformer must be known to check the `incomplete`
        # attribute.
        conformer = polymer.mol.GetConformer(-1).GetId()

        optimizers = [xtb_crude, xtb_normal, xtb_vtight]
        for optimizer in optimizers:
            optimizer.optimize(polymer, conformer)
            if (polymer, conformer) not in in optimizer.incomplete:
                break

    See Also
    --------
    #. https://xtb-docs.readthedocs.io/en/latest/setup.html

    """

    def __init__(self,
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
                 use_cache=False):
        """
        Initializes a :class:`XTB` instance.

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
            molecule and conformer.

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

        self.xtb_path = xtb_path
        self.gfn_version = str(gfn_version)
        self.output_dir = output_dir
        self.opt_level = opt_level
        self.max_runs = max_runs
        self.calculate_hessian = calculate_hessian
        self.num_cores = str(num_cores)
        self.electronic_temperature = str(electronic_temperature)
        self.solvent = solvent
        self.solvent_grid = solvent_grid
        self.charge = str(charge)
        self.num_unpaired_electrons = str(num_unpaired_electrons)
        self.unlimited_memory = unlimited_memory
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

    def _complete(self, output_file):
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
            if self.calculate_hessian:
                return not self._has_neg_frequencies(output_file)
            else:
                return True
        elif os.path.exists('NOT_CONVERGED'):
            raise XTBConvergenceError('Optimization not converged.')
        else:
            raise XTBOptimizerError('Optimization failed to complete')

    def _run_xtb(self, xyz, out_file):
        """
        Runs GFN-xTB.

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
        if self.unlimited_memory:
            memory = 'ulimit -s unlimited ;'
        else:
            memory = ''

        # Set optimization level and type.
        if self.calculate_hessian:
            # Do optimization and check hessian.
            optimization = f'--ohess {self.opt_level}'
        else:
            # Do optimization.
            optimization = f'--opt {self.opt_level}'

        if self.solvent is not None:
            solvent = f'--gbsa {self.solvent} {self.solvent_grid}'
        else:
            solvent = ''

        cmd = (
            f'{memory} {self.xtb_path} {xyz} --gfn {self.gfn_version} '
            f'{optimization} --parallel {self.num_cores} '
            f'--etemp {self.electronic_temperature} '
            f'{solvent} --chrg {self.charge} '
            f'--uhf {self.num_unpaired_electrons}'
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

    def _run_optimizations(self, mol, conformer):
        """
        Run loop of optimizations of the molecule `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`
            The conformer to use.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the calculation is complete and
            ``False`` if the calculation is incomplete.

        """
        for run in range(self.max_runs):
            xyz = f'input_structure_{run+1}.xyz'
            out_file = f'optimization_{run+1}.output'
            mol.write(xyz, conformer=conformer)
            self._run_xtb(xyz=xyz, out_file=out_file)
            # Check if the optimization is complete.
            coord_file = 'xtbhess.coord'
            coord_exists = os.path.exists(coord_file)
            output_xyz = 'xtbopt.xyz'
            opt_complete = self._complete(out_file)
            if not opt_complete:
                if coord_exists:
                    # The calculation is incomplete.
                    # Update mol from xtbhess.coord and continue.
                    mol.update_from_turbomole(
                        path=coord_file,
                        conformer=conformer
                    )
                else:
                    # Update mol from xtbopt.xyz.
                    mol.update_from_xyz(
                        path=output_xyz,
                        conformer=conformer
                    )
                    # If the negative frequencies are small, then GFN
                    # may not produce the restart file. If that is the
                    # case, exit optimization loop and warn.
                    self.incomplete.add((mol, conformer))
                    logging.warning(
                        'Small negative frequencies present in '
                        f'{mol.name} conformer {conformer}.'
                    )
                    return False
            else:
                # Optimization is complete.
                # Update mol from xtbopt.xyz.
                mol.update_from_xyz(
                    path=output_xyz,
                    conformer=conformer
                )
                break
        return opt_complete

    def optimize(self, mol, conformer=-1):
        """
        Optimizes the molecule `mol` using xTB.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        Returns
        -------
        None : :class:`NoneType`

        """

        if conformer == -1:
            conformer = mol.mol.GetConformer(conformer).GetId()
        key = (mol, conformer)

        # Remove (mol, conformer) from self.incomplete if present.
        if key in self.incomplete:
            self.incomplete.remove(key)

        if self.output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self.output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        os.chdir(output_dir)

        try:
            complete = self._run_optimizations(mol, conformer)
        finally:
            os.chdir(init_dir)

        if not complete:
            self.incomplete.add(key)
            logging.warning(
                'Optimization is incomplete for '
                f'{mol.name} conformer {conformer}.'
            )
