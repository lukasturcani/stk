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
    BoundedSelector,
)


logger = logging.getLogger(__name__)


class If(
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
    Selector,
):
    """

    """

    def __init__(
        self,
        condition,
        true_calculator,
        false_calculator,
        **kwargs,
    ):
        """

        """

        self._condition = condition
        self._true_calculator = true_calculator
        self._false_calculator = false_calculator
        super().__init__(**kwargs)

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


class TryCatch(
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
):
    """

    """

    def __init__(
        self,
        try_calculator,
        catch_calculator,
        catch_type=Exception,
        **kwargs,
    ):
        """

        """

        self._try_calculator = try_calculator
        self._catch_calculator = catch_calculator
        self._catch_type = catch_type
        super().__init__(**kwargs)

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

    def _normalize(self, population):
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
    Optimizer,
    FitnessNormalizer,
    BoundedSelector,
    # Putting Selector here is redundant but it makes it clear to
    # users that Sequence is compatible with all Selectors.
    Selector,
):
    """

    """

    def __init__(self, *calculators, **kwargs):
        """

        """

        self._calculators = calculators
        super().__init__(**kwargs)

    def _optimize(self, mol):
        for calculator in self._calculators:
            calculator.optimize(mol)

    def _normalize(self, population):
        for calculator in self._calculators:
            calculator.normalize(population)

    def select(self, population):
        iterables = (
            calculator.select(population)
            for calculator in self._calculators
        )
        yield from it.islice(it.chain(*iterables), self._num_batches)


class Random(
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
    Selector,
):
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
    Crosser,
    EnergyCalculator,
    FitnessCalculator,
    FitnessNormalizer,
    Mutator,
    Optimizer,
    Selector,
):
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
