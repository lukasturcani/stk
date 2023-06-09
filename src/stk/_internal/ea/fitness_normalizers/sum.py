"""
Sum
===

"""

from .fitness_normalizer import FitnessNormalizer


class Sum(FitnessNormalizer):
    """
    Combines fitness value components into a single fitness value.

    Examples
    --------
    *Combining Fitness Value Components Into a Single Value*

    Sometimes, your :class:`.FitnessCalculator` will return a fitness
    value which is a :class:`tuple` of multiple numbers. Each number
    represents a different property of the molecule, which
    contributes to the final fitness value. This fitness normalizer
    combines these fitness value components into a single number, by
    taking their sum.

    .. testcode:: combining-fitness-value-components

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
                fitness_value=(1, -2, 3),
                normalized=False,
            ),
        )

        # Create the normalizer.
        sum_normalizer = stk.Sum()

        normalized_population = tuple(sum_normalizer.normalize(
            population=population,
        ))
        normalized_record, = normalized_population

        assert normalized_record.get_fitness_value() == 2


    *Selectively Normalizing Fitness Values*

    Sometimes, you only want to normalize some members of a population,
    for example if some do not have an assigned fitness value,
    because the fitness calculation failed for whatever reason.
    You can use the `filter` parameter to exclude records from the
    normalization

    .. testcode:: selectively-normalizing-fitness-values

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
                fitness_value=(1, -2, 3),
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

        normalizer = stk.Sum(
            # Only normalize values which are not None.
            filter=lambda population, record:
                record.get_fitness_value() is not None,
        )
        normalized_population = tuple(normalizer.normalize(population))
        normalized_record1, normalized_record2 = normalized_population
        assert normalized_record1.get_fitness_value() == 2
        assert normalized_record2.get_fitness_value() is None


    """

    def __init__(self, filter=lambda population, record: True):
        """
        Initialize a :class:`.Sum` instance.

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
        for record in population:
            if self._filter(population, record):
                yield record.with_fitness_value(
                    fitness_value=sum(record.get_fitness_value()),
                )
            else:
                yield record
