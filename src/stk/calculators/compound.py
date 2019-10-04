"""
Compound Calculators
====================

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


logger = logging.getLogger(__name__)


class If(
    Optimizer,
    EnergyCalculator,
    FitnessCalculator,
    Mutator,
    Crosser,
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

    def _optimize(self, mol):
        if self._condition(mol):
            return self._true_calculator.optimize(mol)
        return self._false_calculator.optimize(mol)

    def _get_energy(self, mol):
        if self._condition(mol):
            return self._true_calculator.get_energy(mol)
        return self._false_calculator.get_energy(mol)

    def _get_fitness(self, mol):
        if self._condition(mol):
            return self._true_calculator.get_fitness(mol)
        return self._false_calculator.get_fitness(mol)

    def _mutate(self, mol):
        if self._condition(mol):
            return self._true_calculator.mutate(mol)
        return self._false_calculator.mutate(mol)

    def _cross(self, *mols):
        if self._condition(*mols):
            return self._true_calculator.cross(*mols)
        return self._false_calculator.cross(*mols)


class TryCatch(
    Optimizer,
    EnergyCalculator,
    FitnessCalculator,
    Mutator,
    Crosser,
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


class Sequence(Optimizer, FitnessNormalizer, Selector):
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

    def _select(self, population):
        iterables = (
            calculator.select(population)
            for calculator in self._calculators
        )
        yield from it.islice(it.chain(*iterables), self._num_batches)


class Random(Optimizer, EnergyCalculator, Mutator, Crosser):
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
        calculator = self._generator.choice(
            a=self._calculators,
            p=self._probabilities,
        )
        self._log_choice(calculator)
        return calculator.optimize(mol)

    def _get_energy(self, mol):
        calculator = self._generator.choice(
            a=self._calculators,
            p=self._probabilities,
        )
        self._log_choice(calculator)
        return calculator.get_energy(mol)

    def _mutate(self, mol):
        calculator = self._generator.choice(
            a=self._calculators,
            p=self._probabilities,
        )
        self._log_choice(calculator)
        return calculator.mutate(mol)

    def _cross(self, *mols):
        calculator = self._generator.choice(
            a=self._calculators,
            p=self._probabilities,
        )
        self._log_choice(calculator)
        return calculator.cross(*mols)

    def _log_choice(self, calculator):
        logger.info(
            f'Random selected {calculator.__class__.__name__}.'
        )
