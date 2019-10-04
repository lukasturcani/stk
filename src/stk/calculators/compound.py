"""
Compound Calculators
====================

"""

import logging
import numpy as np

from .optimization import Optimizer
from .energy import EnergyCalculator


logger = logging.getLogger(__name__)


class If(Optimizer, EnergyCalculator):
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


class TryCatch(Optimizer, EnergyCalculator):
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
        except Exception:
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


class Sequence(Optimizer):
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


class Random(Optimizer, EnergyCalculator):
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

    def _log_choice(self, calculator):
        logger.info(
            f'Random selected {calculator.__class__.__name__}.'
        )
