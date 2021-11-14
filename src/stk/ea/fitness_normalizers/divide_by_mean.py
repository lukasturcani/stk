"""
Divide By Mean
==============

"""

import logging
from functools import partial

import numpy as np

from .fitness_normalizer import FitnessNormalizer

logger = logging.getLogger(__name__)


class DivideByMean(FitnessNormalizer):
    """
    Divides fitness values by the population mean.

    While this function can be used if the fitness value of each
    :class:`.Molecule` in the population is a single
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

    Examples
    --------
    *Selectively Normalizing Fitness Values*

    Sometimes you do not want to normalize all the values in a
    population together. For example, if a failed fitness value
    calculation resulted in some records having a fitness value of
    ``None``, you would want to ignore these records from the
    normalization

    .. testcode:: selectively-normalizing-fitness-values

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
                fitness_value=(1., 2., 3.),
                normalized=False,
            ),
            # This will have a fitness value of None.
            stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(building_block, ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ),
        )

        mean_scaler = stk.DivideByMean(
            # Only normalize values which are not None.
            filter=lambda population, record:
                record.get_fitness_value() is not None
        )
        # Calling mean_scaler.normalize() will return a new
        # population holding the molecule records with normalized
        # fitness values.
        normalized_population = tuple(mean_scaler.normalize(
            population=population,
        ))
        normalized_record1, normalized_record2 = normalized_population
        assert np.all(np.equal(
            normalized_record1.get_fitness_value(),
            (1, 1, 1),
        ))

    """

    def __init__(self, filter=lambda population, record: True):
        """
        Initialize a :class:`.DivideByMean` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes two parameters, first is a :class:`tuple`
            of :class:`.MoleculeRecord` instances,
            and the second is a :class:`.MoleculeRecord`. The
            :class:`callable` returns ``True`` or ``False``. Only
            molecules which return ``True`` will have fitness values
            normalized. By default, all molecules will have fitness
            values normalized.
            The instance passed to the `population` argument of
            :meth:`.normalize` is passed as the first argument, while
            the second argument will be passed every
            :class:`.MoleculeRecord` in it, one at a time.

        """

        self._filter = filter

    def normalize(self, population):
        filtered = filter(
            partial(self._filter, population),
            population,
        )
        mean = np.mean(
            a=[record.get_fitness_value() for record in filtered],
            axis=0,
        )
        logger.debug(f'Means used: {mean}')

        for record in population:
            if self._filter(population, record):
                yield record.with_fitness_value(
                    fitness_value=np.divide(
                        record.get_fitness_value(),
                        mean,
                    )
                )
            else:
                yield record
