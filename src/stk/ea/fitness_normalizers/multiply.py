"""
Multiply
========

"""

import numpy as np

from .fitness_normalizer import FitnessNormalizer


class Multiply(FitnessNormalizer):
    """
    Multiplies the fitness values by some coefficient.

    Examples
    --------
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
                fitness_value=(1, 1, 1),
                normalized=False,
            ),
        )

        normalizer = stk.Multiply((1, 2, 3))
        normalized_population = tuple(normalizer.normalize(population))
        normalized_record, = normalized_population
        assert np.all(np.equal(
            normalized_record.get_fitness_value(),
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
                fitness_value=(1, 1, 1),
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

        normalizer = stk.Multiply(
            coefficient=(1, 2, 3),
            # Only normalize values which are not None.
            filter=lambda population, record:
                record.get_fitness_value() is not None,
        )

        # Calling normalizer.normalize() will return a new
        # population holding the molecule records with normalized
        # fitness values.
        normalized_population = tuple(normalizer.normalize(
            population=population,
        ))
        normalized_record1, normalized_record2 = normalized_population
        assert np.all(np.equal(
            normalized_record1.get_fitness_value(),
            (1, 2, 3),
        ))
        assert normalized_record2.get_fitness_value() is None

    """

    def __init__(
        self,
        coefficient,
        filter=lambda population, record: True,
    ):
        """
        Initialize a :class:`.Multiply` instance.

        Parameters
        ----------
        coefficient : :class:`float` or \
                :class:`tuple` of :class:`float`
            The coefficients each fitness value is multiplied by. Can
            be a single number or multiple numbers, depending on the
            form of the fitness value.

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

        self._coefficient = coefficient
        self._filter = filter

    def normalize(self, population):
        for record in population:
            if self._filter(population, record):
                # Use np.multiply here so that both lists and arrays
                # work properly. If * is used [8]*3 will wrongly make
                # [8, 8, 8].
                yield record.with_fitness_value(
                    fitness_value=np.multiply(
                        self._coefficient,
                        record.get_fitness_value(),
                    ),
                )
            else:
                yield record
