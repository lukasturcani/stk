"""
Power
=====

"""

import numpy as np

from .fitness_normalizer import FitnessNormalizer


class Power(FitnessNormalizer):
    """
    Raises fitness values to some power.

    Examples
    --------
    *Raising Fitness Values to a Power*

    Sometimes you might calculate a property for a molecule, where
    that property indicates a low fitness value. You can use
    :class:`.Power` to raise it to the power of -1 to get your
    final fitness value

    .. testcode:: raising-fitness-values-to-a-power

        import stk

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
                fitness_value=1,
                normalized=False,
            ),
            stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(building_block, ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ).with_fitness_value(
                fitness_value=2,
                normalized=False,
            ),
        )

        normalizer = stk.Power(-1)
        normalized_population = tuple(normalizer.normalize(population))
        normalized_record1, normalized_record2 = normalized_population
        assert normalized_record1.get_fitness_value() == 1
        assert normalized_record2.get_fitness_value() == 0.5

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
                fitness_value=(2, 2, 2),
                normalized=False,
            ),
        )

        normalizer = stk.Power((1, -1, 2))
        normalized_population = tuple(normalizer.normalize(population))
        normalized_record, = normalized_population
        assert np.all(np.equal(
            normalized_record.get_fitness_value(),
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
                fitness_value=(2, 2, 2),
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

        normalizer = stk.Power(
            power=(1, -1, 2),
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
            (2, 0.5, 4),
        ))
        assert normalized_record2.get_fitness_value() is None

    """

    def __init__(
        self,
        power,
        filter=lambda population, record: True,
    ):
        """
        Initialize a :class:`.Power` instance.

        Parameters
        ----------
        power : :class:`float` or \
                :class:`tuple` of :class:`float`
            The power each fitness value is raised to. Can
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

        self._power = power
        self._filter = filter

    def normalize(self, population):
        for record in population:
            if self._filter(population, record):
                yield record.with_fitness_value(
                    fitness_value=np.float_power(
                        record.get_fitness_value(),
                        self._power,
                    ),
                )
            else:
                yield record
