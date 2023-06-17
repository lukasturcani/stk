"""
Null Fitness Normalizer
=======================

"""

from .fitness_normalizer import FitnessNormalizer


class NullFitnessNormalizer(FitnessNormalizer):
    """
    Does nothing.

    This normalizer just yields the molecule records passed to it,
    without changing them in any way.

    """

    def normalize(self, population):
        yield from population
