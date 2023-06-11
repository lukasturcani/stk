"""
Null Optimizer
==============

"""

from .optimizer import Optimizer


class NullOptimizer(Optimizer):
    """
    Does not perform an optimization.

    """

    def optimize(self, state):
        return state
