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
from os.path import join
import shutil
from ...utilities import valid_XTB_solvent, XTBExts

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


class XTBOptimizerFailedError(Exception):
    ...


class XTB(Optimizer):
    """
    Uses GFN-xTB to optimize molecules.

    Notes
    -----
    When running :meth:`optimize`, this calculator changes the
    present working directory with :func:`os.chdir`. The original
    working directory will be restored even if an error is raised so
    unless multi-threading is being used this implementation detail
    should not matter.

    If multi-threading is being used an error could occur if two
    different threads need to know about the current working directory
    as this :class:`.Optimizer` can change it from under them.

    Note that this does not have any impact on multi-processing,
    which should always be safe.

    Furthermore, the :meth:`optimize` calculator will check that the
    structure is adequately optimized by checking for negative frequencies
    after a Hessian calculation. A loop over at the given opt_level will
    be performed by default to obtain an optimized structure. However,
    in the examples we outline how to iterate over opt_levels to increase
    convergence criteria and hopefully obtain an optimized structure.

    Documentation for GFN-xTB available:
    https://xtb-docs.readthedocs.io/en/latest/setup.html

    Attributes
    ----------
    gfnxtb_path : :class:`str`
        The path to the GFN-xTB executable.

    gfn_version : :class:`str`
        Parameterization of GFN-xTB to use.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/basics.html

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    max_count : :class:`int`, optional
        Number of optimizations to attempt in a row to remove negative
        frequencies.

    opt_level : :class:`str`, optional
        Optimization level to use.
        Options:
            crude, sloppy, loose, lax, normal, tight, vtight, extreme
        Definitions of levels:
            https://xtb-docs.readthedocs.io/en/latest/optimization.html

    output_dir : :class:`str`, optional
        The name of the directory into which files generated during
        the optimization are written, if ``None`` then
        :func:`uuid.uuid4` is used.

    num_cores : :class:`int`
        The number of cores for GFN-xTB to use. Requires appropriate setup
        of GFN-xTB by user.

    use_cache : :class:`bool`, optional
        If ``True`` :meth:`optimize` will not run twice on the same
        molecule and conformer.

    mem_ulimit : :class: `bool`, optional
        If ``True`` :meth:`optimize` will be run without constraints on
        the stacksize. If memory issues are encountered, this should be ``True``,
        however this may raise issues on clusters.

    etemp : :class:`int`, optional
        Electronic temperature to use (in K). Defaults to 300K.

    solvent : :class:`str`, optional
        Solvent to use in GBSA implicit solvation method.
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    solvent_grid : :class:`str`, optional
        Grid level to use in SASA calculations for GBSA implicit solvent.
        Options:
            normal, tight, verytight, extreme
        For details:
            https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    multiplicity : :class:`str`, optional
        Number of unpaired electrons.

    charge : :class:`str`, optional
        Formal molecular charge. `-` should be used to indicate sign.

    strict : :class:`bool`, optional
        Whether to use the `--strict` during optimization, which turns all
        internal GFN-xTB warnings into errors.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = GFNXTB('/opt/gfnxtb/xtb')
        gfnxtb.optimize(mol)

    Note that for :class:`.MacroMolecule` objects assembled by ``stk``
    :class:`GFNXTB` should usually be used in a
    :class:`OptimizerSequence`. This is because GFN-xTB only uses
    xyz coordinates as input and so will not recognize the long bonds
    created during assembly. An optimizer which can minimize
    these bonds should be used before :class:`GFNXTB`.

    .. code-block:: python

        bb1 = StructUnit2.smiles_init('NCCNCCN', ['amine'])
        bb2 = StructUnit2.smiles_init('O=CCCC=O', ['aldehyde'])
        polymer = Polymer([bb1, bb2], Linear("AB", [0, 0], 3))

        gfnxtb = OptimizerSequence(
            UFF(),
            GFNXTB(gfnxtb_path='/opt/gfnxtb/xtb',
                   mem_ulimit=True)
        )
        gfnxtb.optimize(polymer)

    Note that energies and other properties of optimized structures can be
    extracted using :class:`GFNXTBEnergy`. For example vibrational frequencies,
    the HOMO-LUMO gap and thermodynamic properties (such as the Total Free
    Energy) can be calculated after an optimization with very tight contraints
    with an implicit solvent (THF). Very tight criteria are required to ensure
    that no negative vibrational frequencies are present.

    .. code-block:: python

        gfnxtb = OptimizerSequence(
            UFF(),
            GFNXTB(gfnxtb_path='/opt/gfnxtb/xtb',
                   mem_ulimit=True,
                   opt_level='verytight',
                   solvent='THF')
        )
        gfnxtb.optimize(polymer)

        gfnxtbenergy = GFNXTBEnergy(gfnxtb_path='/opt/gfnxtb/xtb',
                                    mem_ulimit=True,
                                    free=True,
                                    solvent='THF')
        polymer_properties = gfnxtbenergy.energy(polymer)
        polymer_free_energy = polymer_properties['totalfreeenergy']
        polymer_freq = polymer_properties['frequencies']
        polymer_gap = polymer_properties['HLGap']
    """

    def __init__(self,
                 gfnxtb_path,
                 gfn_version='2',
                 output_dir=None,
                 opt_level='normal',
                 max_count=2,
                 num_cores=1,
                 etemp=300,
                 solvent=None,
                 solvent_grid='normal',
                 charge=None,
                 multiplicity=None,
                 use_cache=False,
                 mem_ulimit=False,
                 strict=True):
        """
        Initializes a :class:`GFNXTB` instance.

        Parameters
        ----------
        gfnxtb_path : :class:`str`
            The path to the GFN-xTB or GFN2-xTB executable.

        gfn_version : :class:`str`
            Parameterization of GFN-xTB to use.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/basics.html

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        opt_level : :class:`str`, optional
            Optimization level to use.
            Options:
                crude, sloppy, loose, lax, normal, tight, vtight, extreme
            Definitions of levels:
                https://xtb-docs.readthedocs.io/en/latest/optimization.html

        max_count : :class:`int`, optional
            Number of optimizations to attempt in a row to remove negative
            frequencies.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        num_cores : :class:`int`
            The number of cores for GFN-xTB to use. Requires appropriate setup
            of GFN-xTB by user.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        mem_ulimit : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without constraints on
            the stacksize. If memory issues are encountered, this should be ``True``,
            however this may raise issues on clusters.

        etemp : :class:`int`, optional
            Electronic temperature to use (in K). Defaults to 300K.

        solvent : :class:`str`, optional
            Solvent to use in GBSA implicit solvation method.
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        solvent_grid : :class:`str`, optional
            Grid level to use in SASA calculations for GBSA implicit solvent.
            Options:
                normal, tight, verytight, extreme
            For details:
                https://xtb-docs.readthedocs.io/en/latest/gbsa.html

        multiplicity : :class:`str`, optional
            Number of unpaired electrons.

        charge : :class:`str`, optional
            Formal molecular charge. `-` should be used to indicate sign.

        strict : :class:`bool`, optional
            Whether to use the `--strict` during optimization, which turns all
            internal GFN-xTB warnings into errors.

        """

        self.gfnxtb_path = gfnxtb_path
        self.gfn_version = gfn_version
        self.output_dir = output_dir
        self.opt_level = opt_level
        self.max_count = max_count
        self.num_cores = str(num_cores)
        self.etemp = str(etemp)
        self.solvent = solvent
        if self.solvent is not None:
            self.solvent = solvent.lower()
            valid_XTB_solvent(gfn_version=self.gfn_version,
                                 solvent=self.solvent)
        self.solvent_grid = solvent_grid
        self.charge = charge
        self.multiplicity = multiplicity
        self.mem_ulimit = mem_ulimit
        self.strict = strict
        self.NOT_OPTIMIZED = False
        super().__init__(use_cache=use_cache)

    def __check_neg_frequencies(self, output_file):
        """
        Check for negative frequencies in output_string.

        Formatting based on latest version of GFN-xTB (190418).
        Example line:
        "eigval :       -0.00    -0.00    -0.00     0.00     0.00     0.00"

        Returns
        -------
        :class:`bool`
            Returns `True` if a negative frequency is present
        """
        # get output file in string
        neg_freq = False
        XTBExt = XTBExts(output_file=output_file)
        value = XTBExt.ext_frequencies()
        # checks for one negative frequency
        if min(value) < 0:
            neg_freq = True
        return neg_freq

    def __check_incomplete(self, output_file):
        """
        Check if GFN-xTB optimization was converged.

        Returns
        -------
        :class:`bool`
            Returns `True` if a negative frequency is present. Raises errors if
            optimization did not converge.
        """
        if output_file is None:
            # no simulation run yet
            return True
        # check for convergence
        if os.path.isfile('.xtboptok'):
            # optimization has converged - however:
            # check for negative frequencies -- if they exist, return True
            # to allow for iteration
            return self.check_neg_frequencies(output_file=output_file)
        # else return errors
        elif os.path.isfile('NOT_CONVERGED'):
            raise XTBOptimizerFailedError(f'Optimization not converged.')
        else:
            raise XTBOptimizerFailedError(f'Optimization failed to complete')

    def __write_and_run_command(self, mol, conformer, count):
        """
        Writes and runs the command for GFN.

        Returns
        -------
        out_file : :class:`str`
            Returns output file of the optimization run by the command.
        """
        # when in output_dir -- use relative paths as GFN does not handle
        # full path in command
        xyz = f'input_structure_{count}.xyz'
        out_file = f'optimization_{count}.output'
        mol.write(xyz, conformer=conformer)
        # modify memory limit
        if self.mem_ulimit:
            cmd = ['ulimit -s unlimited ;']
            # allow multiple shell commands to be run in one subprocess
            shell = True
        else:
            cmd = []
            shell = False
        cmd.append(self.gfnxtb_path)
        cmd.append(xyz)
        # set GFN Parameterization
        if self.gfn_version != '2':
            cmd.append('--gfn')
            cmd.append(self.gfn_version)
        # set optimization level and type
        cmd.append('--ohess')
        cmd.append(self.opt_level)
        # set number of cores
        cmd.append('--parallel')
        cmd.append(self.num_cores)
        # add eletronic temp term
        if self.etemp != '300':
            cmd.append('--etemp')
            cmd.append(self.etemp)
        # write solvent section of cmd
        if self.solvent is not None:
            cmd.append('--gbsa')
            cmd.append(self.solvent)
            if self.solvent_grid != 'normal':
                cmd.append(self.solvent_grid)
        # write charge section of cmd
        if self.charge is not None:
            cmd.append('--chrg')
            cmd.append(self.charge)
        # write multiplicity section of cmd
        if self.multiplicity is not None:
            cmd.append('--uhf')
            cmd.append(self.multiplicity)
        # add strict term
        if self.strict is True:
            cmd.append('--strict')
        cmd = ' '.join(cmd)
        f = open(out_file, 'w')
        # uses the shell if mem_ulimit = True and waits until
        # subproces is complete. This is required to run the mem_ulimit_cmd
        # and GFN calculation in one command, which is then closed, which
        # minimizes the risk of unrestricting the memory limits.
        sp.call(cmd, stdin=sp.PIPE, stdout=f, stderr=sp.PIPE,
                 shell=shell)
        f.close()
        return out_file

    def optimize(self, mol, conformer=-1):
        """
        Optimizes the molecule `mol` using GFN-xTB.

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

        if self.output_dir is None:
            output_dir = str(uuid.uuid4().int)
        else:
            output_dir = self.output_dir
        output_dir = os.path.abspath(output_dir)

        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        os.mkdir(output_dir)
        init_dir = os.getcwd()
        try:
            os.chdir(output_dir)
            run_count = 0
            out_file = None
            while self.__check_incomplete(output_file=out_file):
                run_count += 1
                out_file = self.__write_and_run_command(mol=mol,
                                                        conformer=conformer,
                                                        count=run_count)
                # automated restart from GFN produces xtbopt.coord
                if run_count == 1:
                    output_xyz = join(output_dir, 'xtbopt.xyz')
                    mol.update_from_xyz(path=output_xyz, conformer=conformer)
                else:
                    output_coord = join(output_dir, 'xtbopt.coord')
                    mol.update_from_coord(path=output_coord, conformer=conformer)
                    if run_count == self.max_count:
                        self.NOT_OPTIMIZED = True
                        logging.warning(
                            f'Negative frequencies not removed in max_count optimizations')
                        break
        finally:
            os.chdir(init_dir)
