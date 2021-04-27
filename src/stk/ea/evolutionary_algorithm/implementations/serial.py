"""
Serial Evolutionary Algorithm
=============================

"""

from .implementation import Implementation


class Serial(Implementation):
    """
    A serial implementation of the default evolutionary algorithm.

    """

    def get_generations(self, num_generations):
        yield from self._get_generations(num_generations, map)
