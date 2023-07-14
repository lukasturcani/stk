import typing
from collections.abc import Callable, Iterable
from typing import Any

import numpy as np

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class Multiply(FitnessNormalizer[T]):
    """
    Multiplies the fitness values by some coefficient.

    Examples:

        *Multiplying Fitness Values by a Set of Coefficients*

        In this example, assume that each fitness value consists of a
        :class:`tuple` of numbers, each representing a different property
        of the molecule, and each contributing to the final fitness value.
        The properties can be anything, such as  energy, number of atoms
        or diameter.

        If your final fitness value depends on the combination of these
        properties, you will probably first want to scale them with
        :class:`.DivideByMean`. Once this is done, you may want to
        multiply each property by some coefficient, which reflects its
        relative importance to the final fitness value. For example
        if you multiply the value of one property by ``1`` and another
        by ``2``, the second will contribute twice as much to the
        final fitness value, assuming that you get the final fitness value
        by using the :class:`.Sum` normalizer after :class:`.Multiply`.

        Giving a concrete example

        .. testcode:: multiplying-fitness-values-by-a-set-of-coefficients

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
                record: (1, 1, 1),
            }
            normalizer = stk.Multiply((1, 2, 3))
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert np.all(np.equal(
                normalized_fitness_values[record],
                (1, 2, 3),
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
                record1: (1, 1, 1),
                record2: None,
            }
            normalizer = stk.Multiply(
                coefficient=(1, 2, 3),
                # Only normalize values which are not None.
                filter=lambda fitness_values, record:
                    fitness_values[record] is not None,
            )
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert np.all(np.equal(
                normalized_fitness_values[record1],
                (1, 2, 3),
            ))
            assert normalized_fitness_values[record2] is None
    """

    def __init__(
        self,
        coefficient: float | Iterable[float],
        filter: Callable[
            [dict[T, Any], T], bool
        ] = lambda fitness_values, record: True,
    ) -> None:
        """
        Parameters:
            coefficient (float | list[float]):
                The coefficients each fitness value is multiplied by. Can
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
        if not isinstance(coefficient, int | float):
            coefficient = tuple(coefficient)
        self._coefficient = coefficient
        self._filter = filter

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        return {
            record: np.multiply(self._coefficient, fitness_value)
            if self._filter(fitness_values, record)
            else fitness_value
            for record, fitness_value in fitness_values.items()
        }
