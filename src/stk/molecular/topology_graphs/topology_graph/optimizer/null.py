from ..construction_state import ConstructionState
from .optimizer import Optimizer


class NullOptimizer(Optimizer):
    """
    Does not perform an optimization.

    """

    def optimize(self, state: ConstructionState) -> ConstructionState:
        return state
