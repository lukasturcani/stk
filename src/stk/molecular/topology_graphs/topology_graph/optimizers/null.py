"""
Null Optimizer
==============

"""

from .optimizer import Optimizer
from ....construction_state import ConstructionState


class NullOptimizer(Optimizer):
    """
    Does not perform an optimization.

    """

    def optimize(
        self,
        state: ConstructionState,
    ) -> ConstructionState:
        return state
