import typing
from collections.abc import Callable
from typing import Any

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class Sum(FitnessNormalizer[T]):
    """
    Combines fitness value components into a single fitness value.

    Examples:

        *Combining Fitness Value Components Into a Single Value*

        Sometimes, your :class:`.FitnessCalculator` will return a fitness
        value which is a :class:`tuple` of multiple numbers. Each number
        represents a different property of the molecule, which
        contributes to the final fitness value. This fitness normalizer
        combines these fitness value components into a single number, by
        taking their sum.

        .. testcode:: combining-fitness-value-components

            import stk

            building_block = stk.BuildingBlock('BrCCBr', stk.BromoFactory())
            record = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[building_block],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            fitness_values = {
                record: (1, -2, 3),
            }
            sum_normalizer = stk.Sum()
            normalized_fitness_values = sum_normalizer.normalize(
                fitness_values=fitness_values,
            )
            assert normalized_fitness_values[record] == 2

        *Selectively Normalizing Fitness Values*

        Sometimes, you only want to normalize some members of a population,
        for example if some do not have an assigned fitness value,
        because the fitness calculation failed for whatever reason.
        You can use the `filter` parameter to exclude records from the
        normalization

        .. testcode:: selectively-normalizing-fitness-values

            import stk

            building_block = stk.BuildingBlock('BrCCBr', stk.BromoFactory())
            record1 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[building_block],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            record2 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[building_block],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            fitness_values = {
                record1: (1, -2, 3),
                record2: None,
            }
            normalizer = stk.Sum(
                # Only normalize values which are not None.
                filter=lambda fitness_values, record:
                    fitness_values[record] is not None,
            )
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert normalized_fitness_values[record1] == 2
            assert normalized_fitness_values[record2] is None
    """

    def __init__(
        self,
        filter: Callable[
            [dict[T, Any], T], bool
        ] = lambda fitness_values, record: True,
    ) -> None:
        """
        Parameters:
            filter:
                A function wihch returns ``True`` or ``False``. Only
                molecules which return ``True`` will have fitness values
                normalized. By default, all molecules will have fitness
                values normalized.
                The instance passed to the `fitness_values` argument of
                :meth:`.normalize` is passed as the first argument, while
                the second argument will be passed every
                :class:`.MoleculeRecord` in it, one at a time.
        """
        self._filter = filter

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        return {
            record: sum(fitness_value)
            if self._filter(fitness_values, record)
            else fitness_value
            for record, fitness_value in fitness_values.items()
        }
