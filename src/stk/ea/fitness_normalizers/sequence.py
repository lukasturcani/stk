"""
Normalizer Sequence
===================

"""

from .fitness_normalizer import FitnessNormalizer


class NormalizerSequence(FitnessNormalizer):
    def __init__(self, fitness_normalizers):
        self._fitness_normalizers = fitness_normalizers

    def normalize(self, population):
        for normalizer in self._fitness_normalizers:
            population = tuple(normalizer.normalize(population))
        yield from population
