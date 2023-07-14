import logging
import typing
from collections.abc import Callable
from typing import Any

import numpy as np

from .fitness_normalizer import FitnessNormalizer

logger = logging.getLogger(__name__)

T = typing.TypeVar("T")


class DivideByMean(FitnessNormalizer[T]):
    """
    Divides fitness values by the population mean.

    While this function can be used if the fitness value of each
    molecule in the population is a single
    number, it is most useful when the fitness value is a
    :class:`tuple` of numbers. In this case, it is necessary to somehow
    combine the numbers so that a single fitness value is produced.
    For example, take a fitness value which is the vector holding the
    properties ``[energy, diameter, num_atoms]``. For a given molecule
    these numbers may be something like ``[200,000, 12, 140]``. If we
    were to sum these numbers, the energy term would dominate the final
    fitness value. In order to combine these numbers we can divide them
    by the population averages. For example, if the average energy
    of molecules in the population is ``300,000`` the average diameter
    is ``10`` and the average number of atoms is ``70`` then the
    fitness vector would be scaled to ``[0.5, 1.2, 2]``. These
    numbers are now of a similar magnitude and can be summed to give a
    reasonable value. After division , each value represents how
    much better than the population average each property value is.
    In essence we have removed the units from each parameter.

    Examples:

        *Selectively Normalizing Fitness Values*

        Sometimes you do not want to normalize all the values in a
        population together. For example, if a failed fitness value
        calculation resulted in some records having a fitness value of
        ``None``, you would want to ignore these records from the
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
                record1: (1., 2., 3.),
                record2: None,
            }
            mean_scaler = stk.DivideByMean(
                # Only normalize values which are not None.
                filter=lambda fitness_values, record:
                    fitness_values[record] is not None
            )
            normalized_fitness_values = mean_scaler.normalize(
                fitness_values=fitness_values,
            )
            assert np.all(np.equal(
                normalized_fitness_values[record1],
                (1, 1, 1),
            ))
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
        filtered = {
            record: fitness_value
            for record, fitness_value in fitness_values.items()
            if self._filter(fitness_values, record)
        }
        mean = np.mean(
            a=[fitness_values[record] for record in filtered],
            axis=0,
        )
        logger.debug(f"Means used: {mean}")

        return {
            record: np.divide(fitness_value, mean)
            if self._filter(fitness_values, record)
            else fitness_value
            for record, fitness_value in fitness_values.items()
        }
