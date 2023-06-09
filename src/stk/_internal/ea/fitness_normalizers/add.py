"""
Add
===

"""

import numpy as np

from .fitness_normalizer import FitnessNormalizer


class Add(FitnessNormalizer):
    """
    Adds a number to the fitness values.

    Examples
    --------
    *Incrementing Fitness Values by a Set of Values*

    In this example, assume that each fitness value consists of a
    :class:`tuple` of numbers, each representing a different property
    of the molecule, and each contributing to the final fitness value.
    The properties can be anything, such as  energy, number of atoms
    or diameter.

    Often, if these properties indicate a low fitness value, you
    will take their inverse. However, if these properties can have
    a value of 0, and you try to take their inverse you can end up
    dividing by 0, which is bad. To avoid this, you can add a number,
    like 1, to the fitness values before taking their inverse. This
    normalizer allows you to do this.

    Giving a concrete example

    .. testcode:: incrementing-fitness-values-by-a-set-of-values

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
                fitness_value=(0, 0, 0),
                normalized=False,
            ),
        )

        normalizer = stk.Add((1, 2, 3))
        # Calling normalizer.normalize() will return a new
        # population holding the molecule records with normalized
        # fitness values.
        normalized_population = tuple(normalizer.normalize(
            population=population,
        ))
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

    .. testcode:: incrementing-fitness-values-by-a-set-of-values

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
                fitness_value=(0, 0, 0),
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

        normalizer = stk.Add(
            number=(1, 2, 3),
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
        number,
        filter=lambda population, record: True,
    ):
        """
        Initialize a :class:`.Add` instance.

        Parameters
        ----------
        number : :class:`float` or \
                :class:`tuple` of :class:`float`
            The number each fitness value is increased by. Can
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

        self._number = number
        self._filter = filter

    def normalize(self, population):
        for record in population:
            if self._filter(population, record):
                # Use np.add here so that both tuples and arrays
                # work properly.
                yield record.with_fitness_value(
                    fitness_value=np.add(
                        self._number,
                        record.get_fitness_value(),
                    ),
                )
            else:
                yield record
