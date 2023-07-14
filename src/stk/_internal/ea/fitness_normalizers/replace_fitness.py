import typing
from collections.abc import Callable
from typing import Any

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class ReplaceFitness(FitnessNormalizer[T]):
    """
    Replaces fitness values of a certain value with a new value.

    Examples:
        *Giving a Fitness Value to Failed Calculations*

        You want to replace fitness values which are ``None`` with
        half the minimum fitness value in the population. A fitness value
        may be ``None`` because the fitness calculation failed for some
        reason.

        .. testcode:: giving-a-fitness-value-to-failed-calculations

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
                record1: 2,
                record2: None,
            }

            def get_minimum_fitness_value(fitness_values):
                return min(
                    value for value in fitness_values.values()
                    if value is not None
                )

            replacer = stk.ReplaceFitness(
                # The replacement is half the minimum fitness value in the
                # population.
                get_replacement=lambda fitness_values:
                    get_minimum_fitness_value(fitness_values) / 2,
                # Only replace fitness values which are None.
                filter=lambda fitness_values, record:
                    fitness_values[record] is None,
            )

            # Use the replacer.
            normalized_fitness_values = replacer.normalize(fitness_values)
            assert normalized_fitness_values[record1] == 2
            assert normalized_fitness_values[record2] == 1
    """

    def __init__(
        self,
        get_replacement: Callable[[dict[T, Any]], Any],
        filter: Callable[
            [dict[T, Any], T], bool
        ] = lambda fitness_values, record: True,
    ) -> None:
        """
        Parameters:

            get_replacement:
                Calculates the replacement fitness value. The
                entire population passed to :meth:`.normalize`
                is passed to this parameter, regardless
                of what is passed to the `filter` parameter.

            filter:
                A function which returns ``True`` or ``False``. Only
                molecules which return ``True`` will have fitness values
                normalized. By default, all molecules will have fitness
                values normalized.
                The instance passed to the `fitness_values` argument of
                :meth:`.normalize` is passed as the first argument, while
                the second argument will be passed every
                :class:`.MoleculeRecord` in it, one at a time.
        """
        self._get_replacement = get_replacement
        self._filter = filter

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        replacement = self._get_replacement(fitness_values)
        return {
            record: replacement
            if self._filter(fitness_values, record)
            else fitness_value
            for record, fitness_value in fitness_values.items()
        }
