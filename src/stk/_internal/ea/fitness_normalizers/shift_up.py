import typing
from collections.abc import Callable
from functools import partial
from typing import Any

import numpy as np

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class ShiftUp(FitnessNormalizer[T]):
    """
    Shifts negative fitness values to be positive.

    Assume you have a vector-valued fitness value, where each number
    represents a different property of the molecule, for example,
    ``[1, -10, 1]``

    One way to convert the vector-valued fitness value into a
    scalar fitness value is by summing the elements, and the result in
    this case would be ``-8``. Clearly this doesn't work, because the
    resulting fitness value is not a positive number. To fix this,
    the ``-10`` should be shifted to a positive value.

    :class:`.ShiftUp` finds the minimum value of each element in the
    vector-valued fitness value across the entire population, and for
    elements where this minimum value is less than ``0``, shifts up
    the element value for every molecule in the population, so that the
    minimum value in the entire population is ``1``.

    For example, take a population with the vector-valued fitness
    values

    .. code-block:: python

            [1, -5, 5]
            [3, -10, 2]
            [2, 20, 1]

    After normalization the fitness values will be.

    .. code-block:: python

            [1, 6, 5]
            [3, 1, 2]
            [2, 31, 1]

    :class:`.ShiftUp` also works when the fitness value is a
    single value.

    Examples:

        *Ensuring Positive Fitness Values*

        Here you final fitness value is calculated by taking a
        :class:`.Sum` of the different components of the fitness value.
        To ensure that the final sum is positive, each component must
        also be positive.

        .. testcode:: ensuring-positive-fitness-values

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
                record1: (1, -2, 3),
                record2: (4, 5, -6),
            }

            # Create the normalizer.
            shifter = stk.ShiftUp()

            normalized_fitness_values = shifter.normalize(fitness_values)
            assert np.all(np.equal(
                normalized_fitness_values[record1],
                (1, 1, 10),
            ))
            assert np.all(np.equal(
                normalized_fitness_values[record2],
                (4, 8, 1),
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
                record1: (1, -2, 3),
                record2: None,
            }
            normalizer = stk.ShiftUp(
                # Only normalize values which are not None.
                filter=lambda fitness_values, record:
                    fitness_values[record] is not None,
            )
            normalized_fitness_values = normalizer.normalize(fitness_values)
            assert np.all(np.equal(
                normalized_fitness_values[record1],
                (1, 1, 3),
            ))
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
                A function which returns ``True`` or ``False``. Only
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
        filtered = filter(
            partial(self._filter, fitness_values),
            fitness_values,
        )
        # Get all the fitness arrays in a matrix.
        fmat = np.array([fitness_values[record] for record in filtered])

        # Get the minimum value of each element across the population.
        # keepdims ensures that np.min returns a 1-D array, because
        # it will be True if fitness values are scalar and False if
        # they are array-valued.
        mins = np.min(fmat, axis=0, keepdims=len(fmat.shape) == 1)

        # Convert all elements in mins which are not to be shifted to 0
        # and make the shift equal to the minimum value + 1.
        shift = np.zeros(len(mins))
        for i, min_ in enumerate(mins):
            if min_ <= 0:
                shift[i] = 1 - min_

        return {
            record: fitness_value + shift
            if self._filter(fitness_values, record)
            else fitness_value
            for record, fitness_value in fitness_values.items()
        }
