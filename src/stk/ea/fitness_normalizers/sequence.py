"""
Normalizer Sequence
===================

"""

from .fitness_normalizer import FitnessNormalizer


class NormalizerSequence(FitnessNormalizer):
    """
    Applies other normalizers in sequence.

    Examples
    --------
    *Using Multiple Fitness Normalizers*

    You want to apply multiple fitness normalizations in sequence,
    for example, by first using :class:`.DivideByMean`, followed by
    :class:`.Sum`

    .. testcode:: using-multiple-fitness-normalizers

        import stk
        import numpy as np

        building_block = stk.BuildingBlock(
            smiles='BrCCBr',
            functional_groups=[stk.BromoFactory()],
        )

        population = (
            stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(building_block, ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ).with_fitness_value(
                fitness_value=(1, 100, 1000),
                normalized=False,
            ),
        )

        # Create the normalizer.
        normalizer = stk.NormalizerSequence(
            fitness_normalizers=(
                stk.DivideByMean(),
                stk.Sum(),
            ),
        )
        normalized_population = normalizer.normalize(population)
        normalized_record, = normalized_population
        assert normalized_record.get_fitness_value() == 3

    """

    def __init__(self, fitness_normalizers):
        """
        Initialize a :class:`.NormalizerSequence`.

        Parameters
        ----------
        fitness_normalizers : :class:`iterable`
            The :class:`.FitnessNormalizer` instances which should be
            used in sequence.

        """

        self._fitness_normalizers = tuple(fitness_normalizers)

    def normalize(self, population):
        for normalizer in self._fitness_normalizers:
            population = tuple(normalizer.normalize(population))
        yield from population
