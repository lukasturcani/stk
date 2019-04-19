"""
Defines optimizers.

Optimizers are objects used to optimize molecules. Each optimizer is
initialized with some settings and used to optimize a molecule
with :meth:`Optimizer.optimize`.

.. code-block:: python

    import rdkit.Chem.AllChem as rdkit
    mol = StructUnit2.smiles_init('NCCCN', ['amine'])
    mmff = RDKitForceField(rdkit.MMFFOptimizeMolecule)
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
:class:`OptimizerPipeline` may be used for this.

.. code-block:: python

    # Create a new optimizer which chains the previously defined
    # mmff and etkdg optimizers.
    opt_pipeline = OptimizerPipeline(etkdg, mmff)

    # Run each optimizer in sequence.
    opt_pipeline.optimize(polymer)

By default, running :meth:`Optimizer.optimize` twice in a row will
perform an optimization a second time on a molecule. If we want to
skip optimizations on molecules which have already been optimized
we can use the :attr:`skip_optimized` flag.

.. code-block:: python

    skipping_etkdg = RDKitEmbedder(rdkit.ETKDG(), skip_optimized=True)
    # This does nothing because polymer has already been optimized by
    # the other optimizers.
    skipping_etkdg.optimize(polymer)

When :attr:`skip_optimized` is set to ``True`` the optimizer will
check the :attr:`.Molecule.optimized` flag on the molecule. If the
flag is set to ``True``, then no optimization will take place.

.. _`adding optimizers`:

Extending stk: Adding optimizers.
---------------------------------

New optimizers are added by writing them in this module or in a
separate module in this directory. An optimizer is defined as a new
class which inherits :class:`Optimizer`. The new optimizer class must
define a :meth:`~Optimizer.optimize` method. The method must take 2
arguments a mandatory `mol` argument and an optional `conformer`
argument. :meth:`~Optimizer.optimize` will take a molecule and change
its structure however it likes. Beyond this there are no requirements
for optimizers.

"""

import rdkit.Chem.AllChem as rdkit
from functools import wraps
import numpy as np
import logging
import warnings


logger = logging.getLogger(__name__)


def _add_optimized_toggle(optimize):
    """
    Adds toggling of :attr:`.Molecule.optimized`

    Decorates `optimize` so that after running the attribute
    :attr:`.Molecule.optimized` is set to ``True``

    Parameters
    ----------
    optimize : :class:`function`
        A function which is to have :attr:`.Molecule.optimized`
        toggling added to it.

    Returns
    -------
    :class:`function`
        The decorated function.

    """

    @wraps(optimize)
    def inner(self, mol, conformer=-1):
        r = optimize(self, mol, conformer)
        mol.optimized = True
        return r

    return inner


def _add_skipping(optimize):
    """
    Adds skipping to `optimize`.

    Decorates `optimize` so that before running it checks if the
    :attr:`.Molecule.optimized` attribute is ``True``. If it is,
    then the function is not run.

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
        if self.skip_optimized and mol.optimized:
            logger.info(f'Skipping "{mol.name}".')
            return
        return optimize(self, mol, conformer)

    return inner


class Optimizer:
    """
    A base class for optimizers.

    Attributes
    ----------
    skip_optimized : :class:`bool`
        If ``True`` then :meth:`optimize` returns immediately for
        molecules where :attr:`.Molecule.optimized` is``True``.

    """

    def __init__(self, skip_optimized=False):
        """
        Initializes an :class:`Optimizer`.

        Parameters
        ----------
        skip_optimized : :class:`bool`, optional
            If ``True`` then :meth:`optimize` returns immediately for
            molecules where :attr:`.Molecule.optimized` is``True``.

        """

        self.skip_optimized = skip_optimized

    def __init_subclass__(cls, **kwargs):
        cls.optimize = _add_skipping(cls.optimize)
        cls.optimize = _add_optimized_toggle(cls.optimize)
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


class OptimizerPipeline(Optimizer):
    """
    Chains multiple :class:`Optimizer` instances together.

    Attributes
    ----------
    optimizers : :class:`tuple` of :class:`Optimizer`
        A number of optimizers, each of which gets applied to a
        molecule, based on the order in this :class:`tuple`.

    Examples
    --------
    Let's say we want to embed a molecule with ETKDG first and then
    minimize it with the MMFF force field.

    >>> import rdkit.Chem.AllChem as rdkit
    >>> mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
    >>> etkdg = RDKitEmbedder(rdkit.ETKDG())
    >>> mmff = RDKitForceField(rdkit.MMFFOptimizeMolecule)
    >>> optimizer = OptimizerPipeline(etkdg, mmff)
    >>> optimizer.optimize(mol)

    """

    def __init__(self, *optimizers, skip_optimized=False):
        """
        Initializes a :class:`OptimizerPipeline` instance.

        Parameters
        ----------
        *optimizers : :class:`Optimizer`
            A number of optimizers, each of which gets applied to a
            molecule, based on the order given.

        """

        self.optimizers = optimizers
        super().__init__(skip_optimized=skip_optimized)

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

        for optimizer in self.optimizers:
            cls_name = optimizer.__class__.__name__
            logger.info(f'Using {cls_name} on "{mol.name}".')
            optimizer.optimize(mol, conformer)


class CageOptimizerPipeline(OptimizerPipeline):
    """
    Applies optimizations to cages.

    The difference between this class and :class:`.OptimizerPipeline`
    is that after each optimizer in the pipeline is used, the number
    of windows in the cage is determined. If there are fewer windows
    than expected, it means that the cage is collapsed and the pipeline
    is stopped early.

    """

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


class RaisingOptimizerError(Exception):
    ...


class RaisingOptimizer(Optimizer):
    """
    Raises and optimizes at random.

    This optimizer is used for debugging to simulate optimization
    functions which sometimes completely successfully and sometimes
    randomly fail.

    Attributes
    ----------
    fn : :class:`function`
        When the optimizer does not fail, it uses this optimization
        function to optimize molecules.

    fail_chance : :class:`float`
        The probability that the optimizer will raise an error each
        time :meth:`optimize` is used.

    Examples
    --------
    >>> mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
    >>> etkdg = RDKitEmbedder(rdkit.ETKDG())
    >>> partial_raiser = RaisingOptimizer(etkdg.optimize)
    >>> partial_raiser.optimize(mol)

    """

    def __init__(self, fn, fail_chance=0.5, skip_optimized=False):
        """
        Initializes :class:`PartialRaiser`.

        Parameters
        ----------
        fn : :class:`function`
            When the optimizer does not fail, it uses this optimization
            function to optimize molecules.

        fail_chance : :class:`float`, optional
            The probability that the optimizer will raise an error each
            time :meth:`optimize` is used.

        skip_optimized : :class:`bool`, optional
            If ``True`` then :meth:`optimize` returns immediately for
            molecules where :attr:`.Molecule.optimized` is``True``.

        """

        self.fn = fn
        self.fail_chance = fail_chance
        super().__init__(skip_optimized=skip_optimized)

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
        self.fn(mol)


class RDKitForceField(Optimizer):
    """
    Optimizes the structure of molecules using ``rdkit``.

    Attributes
    ----------
    rdkit_fn : :class:`function`
        An rdkit optimization function, for example
        :func:`rdkit.UFFOptimizeMolecule` or
        :func:`rdkit.MMFFOptimizeMolecule`.

    Examples
    --------
    Use MMFF to optimize a molecule.

    >>> import rdkit.Chem.AllChem as rdkit
    >>> mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
    >>> ff = RDKitForceField(rdkit.MMFFOptimizeMolecule)
    >>> ff.optimize(mol)

    """

    def __init__(self, rdkit_fn, skip_optimized=False):
        """
        Initializes a :class:`RDKitForceField` instance.

        Parameters
        ----------
        rdkit_fn : :class:`function`
            An rdkit optimization function, for example
            :func:`rdkit.UFFOptimizeMolecule` or
            :func:`rdkit.MMFFOptimizeMolecule`.

        skip_optimized : :class:`bool`, optional
            If ``True`` then :meth:`optimize` returns immediately for
            molecules where :attr:`.Molecule.optimized` is``True``.

        """

        self.rdkit_fn = rdkit_fn
        super().__init__(skip_optimized=skip_optimized)

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

        # Sanitize then optimize the rdkit molecule.
        rdkit.SanitizeMol(mol.mol)
        self.rdkit_fn(mol.mol, confId=conformer)


class RDKitEmbedder(Optimizer):
    """
    Uses :func:`rdkit.EmbedMolecule` to find an optimized structure.

    Attributes
    ----------
    params : :class:`rdkit.EmbedParameters`
        The parameters used for the embedding.

    Examples
    --------
    Use ETKDG to generate an optimized structure.

    >>> import rdkit.Chem.AllChem as rdkit
    >>> mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
    >>> embedder = RDKitEmbedder(rdkit.ETKDG())
    >>> embedder.optimize(mol)

    References
    ----------
    .. [#] http://pubs.acs.org/doi/pdf/10.1021/acs.jcim.5b00654

    """

    def __init__(self, params, skip_optimized=False):
        """
        Initializes a :class:`RDKitEmbedder` instance.

        Parameters
        ----------
        params : :class:`rdkit.EmbedParameters`
            The parameters used for the embedding.

        skip_optimized : :class:`bool`, optional
            If ``True`` then :meth:`optimize` returns immediately for
            molecules where :attr:`.Molecule.optimized` is``True``.

        """

        self.params = params
        super().__init__(skip_optimized=skip_optimized)

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
