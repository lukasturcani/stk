import typing
from collections.abc import Iterable
from typing import Any

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class NormalizerSequence(FitnessNormalizer[T]):
    """
    Applies other normalizers in sequence.

    Examples:

        *Using Multiple Fitness Normalizers*

        You want to apply multiple fitness normalizations in sequence,
        for example, by first using :class:`.DivideByMean`, followed by
        :class:`.Sum`

        .. testcode:: using-multiple-fitness-normalizers

            import stk
            import numpy as np

            building_block = stk.BuildingBlock('BrCCBr', stk.BromoFactory())
            record = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[building_block],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            fitness_values = {
                record: (1, 100, 1000),
            }
            # Create the normalizer.
            normalizer = stk.NormalizerSequence(
                fitness_normalizers=[
                    stk.DivideByMean(),
                    stk.Sum(),
                ],
            )
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert normalized_fitness_values[record] == 3
    """

    def __init__(
        self,
        fitness_normalizers: Iterable[FitnessNormalizer[T]],
    ) -> None:
        """
        Parameters:
            fitness_normalizers (list[FitnessNormalizer[T]]):
                The :class:`.FitnessNormalizer` instances which should be
                used in sequence.
        """
        self._fitness_normalizers = tuple(fitness_normalizers)

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        for normalizer in self._fitness_normalizers:
            fitness_values = normalizer.normalize(fitness_values)
        return fitness_values
