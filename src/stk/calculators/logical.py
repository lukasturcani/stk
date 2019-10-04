"""
Logical Calculators
===================

"""

import logging

from .optimization import Optimizer
from .energy import EnergyCalculator
from .ea import FitnessCalculator


logger = logging.getLogger(__name__)


class If(Optimizer, EnergyCalculator, FitnessCalculator):
    """

    """

    def __init__(
        self,
        condition,
        true_calculator,
        false_calculator,
        use_cache=False,
    ):
        """

        """

        self._condition = condition
        self._true_calculator = true_calculator
        self._false_calculator = false_calculator
        super().__init__(use_cache=use_cache)

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


class TryCatch(Optimizer, EnergyCalculator, FitnessCalculator):
    """

    """

    def __init__(
        self,
        try_calculator,
        catch_calculator,
        use_cache=False,
    ):
        """

        """

        self._try_calculator = try_calculator
        self._catch_calculator = catch_calculator
        super().__init__(use_cache=use_cache)

    def _optimize(self, mol):
        try:
            return self._try_calculator.optimize(mol)
        except Exception:
            self._log_failure()
            return self._catch_calculator.optimize(mol)

    def _get_energy(self, mol):
        try:
            return self._try_calculator.get_energy(mol)
        except Exception:
            self._log_failure()
            return self._catch_calculator.get_energy(mol)

    def _get_fitness(self, mol):
        try:
            return self._try_calculator.get_fitness(mol)
        except Exception:
            self._log_failure()
            return self._catch_calculator.get_fitness(mol)

    def _log_failure(self):
        try_name = self._try_calculator.__class__.__name__
        catch_name = self._catch_calculator.__class__.__name__
        logger.error(
            f'{try_name} failed, trying {catch_name}.',
            exc_info=True
        )
