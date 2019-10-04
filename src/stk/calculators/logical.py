"""
Logical Calculators
===================

"""


from .optimization import Optimizer
from .energy import EnergyCalculator


class If(Optimizer, EnergyCalculator):
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
