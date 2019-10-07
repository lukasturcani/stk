"""
Compound Calculators
====================

#. :class:`.If`
#. :class:`.TryCatch`
#. :class:`.Sequence`
#. :class:`.Random`
#. :class:`.RaisingCalculator`

Compound calculators are calculators which are initialized with other
calculators. They then control how and when those calculators are used.

Compound calculators are special, because they take on the class
of the calculators they hold. This means that if you initialize a
compound calculator of a certain class it will be indistinguishable
from any other calculator of that class. For example, if you
initialize a :class:`.Sequence` with :class:`.Optimizer` instances,
it will be usable in the same way as any other :class:`.Optimizer`

.. code-block:: python

    import stk

    sequence = stk.Sequence(
        stk.ETKDG(),
        stk.MMFF(),
    )
    mol = stk.BuildingBlock('NCCN')
    sequence.optimize(mol)

Not all compound calculators can be used with every calculator type.
A compound calculator is only compatible with classes with which it
shares a base class. So if a compound calculator does not inherit
:class:`.EnergyCalculator`, it cannot be used with instances of
:class:`.EnergyCalculator`. However, if a compound calculator does
inherit :class:`.EnergyCalculator`, it can be initialized with
instances of :class:`.EnergyCalculator` and therefore it can also be
used wherever an :class:`.EnergyCalculator` is required.

If you really want to build up complexity, compound calculators
can hold other compound calculators. For example, this is an
:class:`Optimizer`, which will try to optimize a molecule with one
:class:`.Sequence` and if it raises an error, it will try to
optimize the molecule with a different :class:`.Sequence`

.. code-block:: python

    optimizer = stk.TryCatch(
        try_calculator=stk.Sequence(stk.ETKDG(), stk.MMFF()),
        catch_calculator=stk.Sequence(stk.ETKDG(), stk.UFF()),
    )
    # Try to optimize with ETKDG followed by MMFF first, and if that
    # raises an error, try to optimize with ETKDG followed by
    # UFF.
    optimizer.optimize(mol)


"""

import logging
import numpy as np

from .optimization import Optimizer
from .energy import EnergyCalculator
from .base_calculators import _MoleculeCalculator


logger = logging.getLogger(__name__)


class If(
    _MoleculeCalculator,
    EnergyCalculator,
    Optimizer,
):
    """
    Use a `condition` to pick a calculator.


    Examples
    --------
    Use ETKDG to optimize the structure of a molecule if it has
    less than 200 atoms, else use the MMFF force field

    .. code-block:: python

        import stk

        optimizer = stk.If(
            condition=lambda mol: len(mol.atoms) < 200,
            true_calculator=stk.ETKDG(),
            false_calculator=stk.MMFF(),
        )

        # ETKDG will be run on this molecule.
        small_mol = stk.BuildingBlock('NCCN', ['amine'])
        optimizer.optimize(small_mol)

        # MMFF will be run on polymer.
        polymer = stk.ConstructedMolecule(
            building_blocks=[bb],
            topology_graph=stk.polymer.Linear('A', 500),
        )
        optimizer.optimize(polymer)

    """

    def __init__(
        self,
        condition,
        true_calculator,
        false_calculator,
        use_cache=False,
    ):
        """
        Initialize a :class:`.If` instance.

        Parameters
        ----------
        condition : :class:`callable`
            A callable which is passed the arguments which are to
            be passed to the calculators. If it returns ``True``
            then `true_calculator` will be used, if ``False``, then
            `false_calculator` will be used.

        true_calculator : see base classes
            The calculator to use if `condition` returns ``True``.

        false_calculator : see base classes
            The calculator to use if `condition` returns ``False``.

        use_cache : :class:`bool`, optional
            When :class:`.If` is used as a
            :class:`.MoleculeCalculator`, this toggles the results
            cache. When :class:`.If` is used as an
            :class:`EAOperation`, this toggles use of the molecular
            cache.

        """

        self._condition = condition
        self._true_calculator = true_calculator
        self._false_calculator = false_calculator
        self._use_cache = use_cache

    def _optimize(self, mol):
        if self._condition(mol):
            return self._true_calculator.optimize(mol)
        return self._false_calculator.optimize(mol)

    def _get_energy(self, mol):
        if self._condition(mol):
            return self._true_calculator.get_energy(mol)
        return self._false_calculator.get_energy(mol)


class TryCatch(
    _MoleculeCalculator,
    EnergyCalculator,
    Optimizer,
):
    """
    Try one calculator and if it raises an error use another.

    Examples
    --------
    Try using the MMFF force field to optimize a molecule and if it
    fails use the UFF force field

    .. code-block:: python

        import stk

        optimizer = stk.TryCatch(
            try_calculator=stk.MMFF(),
            catch_calculator=stk.UFF(),
        )
        mol = stk.BuildingBlock('NCCN')
        optimizer.optimize(mol)

    Try using the MMFF force field to calculate energy and if it
    fails use the UFF force field

    .. code-block:: python

        energy_calculator = stk.TryCatch(
            try_calculator=stk.MMFFEnergy(),
            catch_calculator=stk.UFFEnergy(),
        )
        mol = stk.BuildingBlock('NCCN')
        energy = energy_calculator.get_energy(mol)

    """

    def __init__(
        self,
        try_calculator,
        catch_calculator,
        catch_type=Exception,
        use_cache=True,
    ):
        """
        Initialize a :class:`.TryCatch` instance.

        Parameters
        ----------
        try_calculator : see base classes
            The calculator to try using first.

        catch_calculator : see base classes
            The calculator to use if `try_calculator` raises an
            error.

        catch_type : :class:`Exception`
            The `catch_calculator` will only be run if the
            `try_calculator` raises an exception of this type.

        use_cache : :class:`bool`, optional
            When used as a :class:`.MoleculeCalculator`, this toggles
            use of the results cache. When is used as an
            :class:`EAOperation`, this toggles use of the molecular
            cache.

        """

        self._try_calculator = try_calculator
        self._catch_calculator = catch_calculator
        self._catch_type = catch_type
        self._use_cache = use_cache

    def _optimize(self, mol):
        try:
            return self._try_calculator.optimize(mol)
        except self._catch_type:
            self._log_failure()
            return self._catch_calculator.optimize(mol)

    def _get_energy(self, mol):
        try:
            return self._try_calculator.get_energy(mol)
        except self._catch_type:
            self._log_failure()
            return self._catch_calculator.get_energy(mol)

    def _log_failure(self):
        try_name = self._try_calculator.__class__.__name__
        catch_name = self._catch_calculator.__class__.__name__
        logger.error(
            f'{try_name} failed, trying {catch_name}.',
            exc_info=True
        )


class Sequence(
    Optimizer,
):
    """
    Use calculators in sequence.
    """

    def __init__(self, *calculators, **kwargs):
        """

        """

        self._calculators = calculators
        super().__init__(**kwargs)

    def _optimize(self, mol):
        for calculator in self._calculators:
            calculator.optimize(mol)


class Random(EnergyCalculator, Optimizer):
    def __init__(
        self,
        *calculators,
        probabilities=None,
        random_seed=None,
        **kwargs,
    ):
        self._calculators = calculators
        self._probabilities = probabilities
        self._generator = np.random.RandomState(random_seed)
        super().__init__(**kwargs)

    def _optimize(self, mol):
        return self._get_calculator().optimize(mol)

    def _get_energy(self, mol):
        return self._get_calculator().get_energy(mol)

    def _get_calculator(self):
        calculator = self._generator.choice(
            a=self._calculators,
            p=self._probabilities,
        )
        self._log_choice(calculator)
        return calculator

    def _log_choice(self, calculator):
        logger.info(
            f'Random selected {calculator.__class__.__name__}.'
        )


class RaisingCalculatorError(Exception):
    ...


class RaisingCalculator(Optimizer, EnergyCalculator):
    """

    """

    def __init__(
        self,
        calculator,
        fail_chance=0.5,
        random_seed=None,
        **kwargs
    ):
        """

        """

        self._calculator = calculator
        self._fail_chance = fail_chance
        self._generator = np.random.RandomState(random_seed)
        super().__init__(**kwargs)

    def _try_raising(self):
        if self._generator.rand() < self._fail_chance:
            raise RaisingCalculatorError()

    def _optimize(self, mol):
        self._try_raising()
        return self._calculator.optimize(mol)

    def _get_energy(self, mol):
        self._try_raising()
        return self._calculator.get_energy(mol)
