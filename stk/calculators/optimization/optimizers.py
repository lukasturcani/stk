"""
Defines optimizers.

Optimizers are objects used to optimize molecules. Each optimizer is
initialized with some settings and used to optimize a molecule
with :meth:`~Optimizer.optimize`.

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit
    mol = StructUnit2.smiles_init('NCCCN', ['amine'])
    mmff = MMFF()
    mmff.optimize(mol)

    # Optionally, a conformer can be provided.
    mmff.optimize(mol, conformer=2)

    # Optimizers also work with MacroMolecule objects.
    polymer = Polymer([mol], Linear('A', [0], n=3))
    etkdg = RDKitEmbedder(rdkit.ETKDG())
    etkdg.optimize(polymer)

Sometimes it is desirable to chain multiple optimizations, one after
another. For example, before running an optimization, it may be
desirable to embed a molecule first, to generate an initial structure.
:class:`OptimizerSequence` may be used for this.

.. code-block:: python

    # Create a new optimizer which chains the previously defined
    # mmff and etkdg optimizers.
    optimizer_sequence = OptimizerSequence(etkdg, mmff)

    # Run each optimizer in sequence.
    optimizer_sequence.optimize(polymer)

By default, running :meth:`Optimizer.optimize` twice in a row will
perform an optimization a second time on a molecule. If we want to
skip optimizations on molecules which have already been optimized
we can use the :attr:`use_cache` flag.

.. code-block:: python

    caching_etkdg = RDKitEmbedder(rdkit.ETKDG(), use_cache=True)
    # First optimize call runs an optimization.
    caching_etkdg.optimize(polymer)
    # Second call does nothing.
    caching_etkdg.optimize(polymer)

Caching is done on a per :class:`Optimizer` basis. Just because the
molecule has been cached by one :class:`Optimizer` instance does not
mean that a different :class:`Optimizer` instance will no longer
optimize the molecule.

.. _`adding optimizers`:

Extending stk: Making new optimizers.
-------------------------------------

New optimizers can be made by simply making a class which inherits the
:class:`Optimizer` class. In addition to this, the new class must
define a :meth:`~Optimize.optimize` method. The method must take 2
arguments a mandatory `mol` argument and an optional `conformer`
argument. :meth:`~Optimizer.optimize` will take a molecule and change
its structure however it likes. Beyond this there are no requirements.

"""

import logging
import numpy as np
import rdkit.Chem.AllChem as rdkit
import warnings
from glob import iglob
import os
from functools import wraps
import subprocess as sp
import uuid

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
    collapsed, the optiimzation sequence ends immediately.

    Attributes
    ----------
    optimizers : :class:`tuple` of :class:`Optimizer`
        The :class:`Optimizer` objects which are used to optimize a
        :class:`.Cage` molecule.

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


class RDKitEmbedder(Optimizer):
    """
    Uses :func:`rdkit.EmbedMolecule` to find an optimized structure.

    Attributes
    ----------
    params : :class:`rdkit.EmbedParameters`
        The parameters used for the embedding.

    Examples
    --------
    Use ETKDG [#]_ to generate an optimized structure.

    .. code-block:: python

        import rdkit.Chem.AllChem as rdkit
        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        embedder = RDKitEmbedder(rdkit.ETKDG())
        embedder.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, params, use_cache=False):
        """
        Initializes a :class:`RDKitEmbedder` instance.

        Parameters
        ----------
        params : :class:`rdkit.EmbedParameters`
            The parameters used for the embedding.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        """

        self.params = params
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

        conf_id = rdkit.EmbedMolecule(mol.mol, self.params)
        new_conf = rdkit.Conformer(mol.mol.GetConformer(conf_id))
        mol.mol.RemoveConformer(conf_id)
        mol.mol.RemoveConformer(conformer)
        new_conf.SetId(conformer)
        mol.mol.AddConformer(new_conf)


class GFNXTB(Optimizer):
    """
    Uses GFN-xTB to optimize molecules.

    Attributes
    ----------
    gfnxtb_path : :class:`str`
        The path to the GFN-xTB executable.

    num_cores : :class:`str`
        The number of cores for GFN-xTB to use.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        gfnxtb = GFNXTB('/opt/gfnxtb/xtb')
        gfnxtb.optimize(mol)

    """

    def __init__(self, gfnxtb_path, output_dir=None, free=False, opt_level=1,
                 num_cores=1, solvent=None, charge=None, use_cache=False):
        """
        Initializes a :class:`GFNXTB` instance.

        Parameters
        ----------
        gfnxtb_path : :class:`str`
            The path to the GFN-xTB or GFN2-xTB executable.

        output_dir : :class:`str`, optional
            The name of the directory into which files generated during
            the optimization are written, if ``None`` then
            :func:`uuid.uuid4` is used.

        free : :class:`bool`, optional
            If ``True`` :meth:`optimize` will perform a numerical hessian
            calculation on the optimized structure to give Free energy also.

        opt_level : :class:`int`, optional
            Optimization level to use.
            -2 =  , -1 = , 0 = , 1 = , 2 = .

        num_cores : :class:`int`
            The number of cores for GFN-xTB to use.

        use_cache : :class:`bool`, optional
            If ``True`` :meth:`optimize` will not run twice on the same
            molecule and conformer.

        solvent : :class:`str`, optional
            XXXX

        charge : :class:`str`, optional
            XXXX
        """

        self.gfnxtb_path = gfnxtb_path
        self.free = free
        self.opt_level = opt_level
        self.num_cores = str(num_cores)
        self.solvent = solvent
        if self.solvent is not None:
            self.valid_solvent()
        self.charge = charge
        super().__init__(use_cache=use_cache)

    def valid_solvent(self):
        '''Check if solvent is valid.

        '''
        valid_solvents = ['ch2cl2']
        if self.solvent in valid_solvents:
            return True
        else:
            raise(f'{self.solvent} is an invalid solvent!')

    def move_generated_GFN_files(self, basename, output_dir):
        '''
        Move files generated by GFN to output_dir.

        Files and purpose:
            basename.xyz: input xyz file
            basename.output: output GFN log file
            .xtboptok: empty file flagging successful optimization
            xtbopt.log: structure at each step of optimization in xyz format
            xtbopt.xyz: optimized structure as xyz file
            xtbrestart: restart file for optimization
            ??
            ??
            ??

        '''
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        # iterate over files with basename
        for filename in iglob(f'{basename}*'):
            # Do not move the output_dir.
            if filename == output_dir:
                continue
            os.rename(filename, f'{output_dir}/{filename}')

        # gfn specific files
        list_of_files = ['.xtboptok', ]
        for file in list_of_files:
            if os.path.isfile(file):
                os.rename(file, f'{output_dir}/{file}')

    def check_complete(self, basename, output_dir):
        if os.path.isfile(f'{output_dir}/.xtboptok'):
            return True
        else:
            return False

    def optimize(self, mol, conformer=-1, mem_ulimit=False):
        """
        Optimizes a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be optimized.

        conformer : :class:`int`, optional
            The conformer to use.

        mem_ulimit : :class: `bool`, optional
            If ``True`` :meth:`optimize` will be run without constraints on
            the stacksize. Generally should be ``True``, however this may raise
            issues on clusters.

        Returns
        -------
        self.check_complete() : :class:`bool`
            Output from process. True if optimization is successful.

        """

        basename = uuid.uuid4().int

        if self.output_dir is None:
            output_dir = basename
        else:
            output_dir = self.output_dir

        xyz = f'{basename}.xyz'
        out_file = f'{basename}.output'
        mol.write(xyz, conformer=conformer)

        # write solvent section of cmd
        if self.solvent is None:
            solvent_part = ''
        else:
            solvent_part = '--gbsa ' + self.solvent

        # write charge section of cmd
        if self.charge is None:
            charge_part = ''
        else:
            charge_part = '--chrg ' + self.solvent

        # set optimization level
        if self.free is False:
            opt_level_part = '-opt ' + self.opt_level
        elif self.free is True:
            opt_level_part = '-ohess ' + self.opt_level

        # modify memory limit
        if mem_ulimit:
            mem_ulimt_cmd = 'ulimit -s unlimited'+' ;'
            # allow multiple shell commands to be run in one subprocess
            shell = True
        else:
            mem_ulimt_cmd = ''
            shell = False

        cmd = [
            mem_ulimt_cmd,
            self.gfnxtb_path, xyz, opt_level_part,
            '--parallel', self.num_cores,
            solvent_part, charge_part, f'> {out_file}'
        ]

        proc = sp.Popen(
            cmd,
            stdin=sp.PIPE,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            shell=shell,
        )
        output, err = proc.communicate()
        self.move_generated_GFN_files(basename=basename, output_dir=output_dir)
        return self.check_complete(basename=basename, output_dir=output_dir)
