import typing
from collections.abc import Callable, Iterable
from typing import Any

import numpy as np

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class Power(FitnessNormalizer[T]):
    """
    Raises fitness values to some power.

    Examples:

        *Raising Fitness Values to a Power*

        Sometimes you might calculate a property for a molecule, where
        that property indicates a low fitness value. You can use
        :class:`.Power` to raise it to the power of -1 to get your
        final fitness value

        .. testcode:: raising-fitness-values-to-a-power

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
                record1: 1,
                record2: 2,
            }
            normalizer = stk.Power(-1)
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert normalized_fitness_values[record1] == 1
            assert normalized_fitness_values[record2] == 0.5

        *Raising Fitness Values by a Set of Powers*

        In this example, assume that each fitness value consists of a
        :class:`tuple` of numbers, each representing a different property
        of the molecule, and each contributing to the final fitness value.
        The properties can be anything, such as  energy, number of atoms
        or diameter.

        If your final fitness value depends on the combination of these
        properties, you will probably first want to scale them with
        :class:`.DivideByMean`. Once this is done, you may want to
        raise each property by some power. For example
        if you raise the value of one property by ``1`` and another
        by ``-1``, it means that when the value of property 1 is big,
        the fitness value should also be big, but if the value of property
        2 is big, the fitness value should be small.

        Giving a concrete example

        .. testcode:: raising-fitness-values-by-a-set-of-powers

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
                record: (2, 2, 2),
            }
            normalizer = stk.Power((1, -1, 2))
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert np.all(np.equal(
                normalized_fitness_values[record],
                (2, 0.5, 4),
            ))

        *Selectively Normalizing Fitness Values*

        Sometimes, you only want to normalize some members of a population,
        for example if some do not have an assigned fitness value,
        because the fitness calculation failed for whatever reason.
        You can use the `filter` parameter to exclude records from the
        normalization

        .. testcode:: selectively-normalizing-fitness-values

            import stk
            import numpy as np

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
                record1: (2, 2, 2),
                record2: None,
            }
            normalizer = stk.Power(
                power=(1, -1, 2),
                # Only normalize values which are not None.
                filter=lambda fitness_values, record:
                    fitness_values[record] is not None,
            )
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert np.all(np.equal(
                normalized_fitness_values[record1],
                (2, 0.5, 4),
            ))
            assert normalized_fitness_values[record2] is None
    """

    def __init__(
        self,
        power: float | Iterable[float],
        filter: Callable[
            [dict[T, Any], T], bool
        ] = lambda fitness_values, record: True,
    ) -> None:
        """
        Parameters:
            power (float | list[float]):
                The power each fitness value is raised to. Can
                be a single number or multiple numbers, depending on the
                form of the fitness value.

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
        if not isinstance(power, int | float):
            power = tuple(power)
        self._power = power
        self._filter = filter

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        return {
            record: np.float_power(fitness_value, self._power)
            if self._filter(fitness_values, record)
            else fitness_value
            for record, fitness_value in fitness_values.items()
        }
