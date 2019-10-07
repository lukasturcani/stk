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
import itertools as it

from .optimization import Optimizer
from .energy import EnergyCalculator
from .ea import (
    FitnessCalculator,
    Mutator,
    Crosser,
    FitnessNormalizer,
    Selector,
)
from .base_calculators import _MoleculeCalculator, MoleculeCalculator


logger = logging.getLogger(__name__)


class If(
    _MoleculeCalculator,
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
    Selector,
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
            cache.

        """

        self._condition = condition
        self._true_calculator = true_calculator
        self._false_calculator = false_calculator
        self._use_cache = use_cache

    def _cross(self, *mols):
        if self._condition(*mols):
            return self._true_calculator.cross(*mols)
        return self._false_calculator.cross(*mols)

    def _get_energy(self, mol):
        if self._condition(mol):
            return self._true_calculator.get_energy(mol)
        return self._false_calculator.get_energy(mol)

    def _get_fitness(self, mol):
        if self._condition(mol):
            return self._true_calculator.get_fitness(mol)
        return self._false_calculator.get_fitness(mol)

    def _normalize(self, population):
        if self._condition(population):
            return self._true_calculator.normalize(population)
        return self._false_calculator.normalize(population)

    def _mutate(self, mol):
        if self._condition(mol):
            return self._true_calculator.mutate(mol)
        return self._false_calculator.mutate(mol)

    def _optimize(self, mol):
        if self._condition(mol):
            return self._true_calculator.optimize(mol)
        return self._false_calculator.optimize(mol)

    def select(self, population):
        if self._condition(population):
            return self._true_calculator.select(population)
        return self._false_calculator.select(population)

    def set_cache_use(self, use_cache):
        # When the contained calculators are MoleculeCalculators,
        # use the implementation provided by _MoleculeCalculator.
        # This will mean that the use_cache option on the contained
        # calculators will not be changed. Instead If will use its own
        # cache to return previously calculated results.
        if isinstance(self._true_calculator, MoleculeCalculator):
            return super().set_cache_use(use_cache)

        # If the contained calculators are EAOperations then set their
        # use_cache option to what is passed here. This is because
        # the If instance does not control the molecular cache.
        # So the contained calculators need to be made aware of this
        # change directly.
        self._true_calculator.set_cache_use(use_cache)
        self._false_calculator.set_cache_use(use_cache)
        return self


class TryCatch(
    _MoleculeCalculator,
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
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
            use of the results cache.

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

    def _get_fitness(self, mol):
        try:
            return self._try_calculator.get_fitness(mol)
        except self._catch_type:
            self._log_failure()
            return self._catch_calculator.get_fitness(mol)

    def _normalize(self, population, fitness_values):
        try:
            return self._try_calculator.normalize(population)
        except self._catch_type:
            self._log_failure()
            return self._catch_calculator.normalize(population)

    def _mutate(self, mol):
        try:
            return self._try_calculator.mutate(mol)
        except self._catch_type:
            self._log_failure()
            return self._catch_calculator.mutate(mol)

    def _cross(self, *mols):
        try:
            return self._try_calculator.cross(*mols)
        except self._catch_type:
            self._log_failure()
            return self._catch_calculator.cross(*mols)

    def _log_failure(self):
        try_name = self._try_calculator.__class__.__name__
        catch_name = self._catch_calculator.__class__.__name__
        logger.error(
            f'{try_name} failed, trying {catch_name}.',
            exc_info=True
        )


class Sequence(
    _MoleculeCalculator,
    Optimizer,
    FitnessNormalizer,
    Selector,
):
    """
    Use calculators in sequence.

    Examples
    --------
    First use :class:`.ETKDG` to optimize a molecule and then use
    :class:`.UFF`

    .. code-block:: python

        import stk

        optimizer = stk.Sequence(
            stk.ETKDG(),
            stk.UFF(),
        )

        mol = stk.BuildingBlock('NCCN')
        optimizer.optimize(mol)


    """

    def __init__(self, *calculators, use_cache=False):
        """
        Initialize a :class:`.Sequence` instance.

        Parameters
        ----------
        calculators : see base classes
            The calculators to be applied in sequence.


        use_cache : :class:`bool`, optional
            When used as a :class:`.MoleculeCalculator`, this toggles
            use of the results cache.

        """

        self._calculators = calculators
        self._use_cache = use_cache

    def _optimize(self, mol):
        for calculator in self._calculators:
            calculator.optimize(mol)

    def _normalize(self, population, fitness_values):
        # Create a list as population is an iterable which gets used
        # multiple times.
        population = list(population)
        normalized = dict(fitness_values)
        for calculator in self._calculators:
            normalized = calculator.normalize(population, normalized)
        yield from normalized.items()

    def select(self, population):
        iterables = (
            calculator.select(population)
            for calculator in self._calculators
        )
        yield from it.islice(it.chain(*iterables), self._num_batches)


class Random(
    _MoleculeCalculator,
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
    Selector,
):
    """
    Pick a calculator to use at random.

    Examples
    --------
    Pick a random mutation operation

    .. code-block:: python

        mutator = stk.Random(
            stk.RandomBuildingBlock(...),
            stk.RandomTopologyGraph(...),
        )
        mol = stk.ConstructedMolecule(...)

        # Mutate mol with either RandomBuildingBlock or
        # RandomTopologyGraph, at random.
        mutant = mutator.mutate(mol)

    """

    def __init__(
        self,
        *calculators,
        probabilities=None,
        random_seed=None,
        use_cache=False,
    ):
        """
        Initialize a :class:`.Random` instance.

        Parameters
        ----------
        calculators : see base classes
            The calculator, one of which is picked at random each
            time a calculation is requested.

        probabilities : :class:`tuple` of :class:`float`, optional
            The probability of picking each calculator in
            `calculators`. If ``None``, all calculators have an
            equal chance of being picked.

        random_seed : :class:`int`, optional
            The random seed for picking the calculator.

        use_cache : :class:`bool`, optional
            When used as a :class:`.MoleculeCalculator`, this toggles
            use of the results cache.

        """

        self._calculators = calculators
        self._probabilities = probabilities
        self._generator = np.random.RandomState(random_seed)
        self._use_cache = use_cache

    def _optimize(self, mol):
        return self._get_calculator().optimize(mol)

    def _get_energy(self, mol):
        return self._get_calculator().get_energy(mol)

    def _get_fitness(self, mol):
        return self._get_calculator().get_fitness(mol)

    def _normalize(self, population):
        return self._get_calculator().normalize(population)

    def _mutate(self, mol):
        return self._get_calculator().mutate(mol)

    def _cross(self, *mols):
        return self._get_calculator().cross(*mols)

    def select(self, population):
        return self._get_calculator().select(population)

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


class RaisingCalculator(
    _MoleculeCalculator,
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
    Selector,
):
    """
    Raise an error at random or use another calculator.

    This calculator is mainly used for debugging.

    .. code-block:: python

        import stk

        energy_calculator = stk.RaisingCalculator(stk.UFFEnergy())
        mol = stk.BuildingBlock('NCCN')
        # 50% chance of getting the energy with UFF and 50% chance
        # of raising an error.
        energy = energy_calculator.get_energy(mol)

    """

    def __init__(
        self,
        calculator,
        fail_chance=0.5,
        random_seed=None,
        use_cache=False,
    ):
        """
        Initialize a :class:`.RaisingCalculator` instance.

        Parameters
        ----------
        calculator : see base classes
            The calculator to use when an error is not raised.

        fail_chance : :class:`float`, optional
            The probability of raising an error.

        random_seed : :class:`int`, optional
            The random seed to use for raising an error.

        use_cache : :class:`bool`, optional
            When used as a :class:`.MoleculeCalculator`, this toggles
            use of the results cache.

        """

        self._calculator = calculator
        self._fail_chance = fail_chance
        self._generator = np.random.RandomState(random_seed)
        self._use_cache = use_cache

    def _try_raising(self):
        if self._generator.rand() < self._fail_chance:
            raise RaisingCalculatorError()

    def _optimize(self, mol):
        self._try_raising()
        return self._calculator.optimize(mol)

    def _get_energy(self, mol):
        self._try_raising()
        return self._calculator.get_energy(mol)

    def _cross(self, *mols):
        self._try_raising()
        return self._calculator.cross(*mols)

    def _get_fitness(self, mol):
        self._try_raising()
        return self._calculator.get_fitness(mol)

    def _normalize(self, population):
        self._try_raising()
        return self._calculator.normalize(population)

    def _mutate(self, mol):
        self._try_raising()
        return self._calculator.mutate(mol)

    def select(self, population):
        self._try_raising()
        return self._calculator.select(population)
