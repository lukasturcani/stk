"""
Replace Fitness
===============

"""

from .fitness_normalizer import FitnessNormalizer


class ReplaceFitness(FitnessNormalizer):
    """
    Replaces fitness values of a certain value with a new value.

    Examples
    --------
    *Giving a Fitness Value to Failed Calculations*

    You want to replace fitness values which are ``None`` with
    half the minimum fitness value in the population. A fitness value
    may be ``None`` because the fitness calculation failed for some
    reason.

    .. code-block:: python

        import stk

        def get_minimum_fitness_value(population):
            return min(
                record.get_fitness_value() for record in population
                if record.get_fitness_value() is not None
            )

        replacer = stk.ReplaceFitness(
            # The replacement is half the minimum fitness value in the
            # population.
            get_replacement=lambda population:
                get_minimum_fitness_value(population) / 2,
            # Only replace fitness values which are None.
            filter=lambda population, record:
                record.get_fitness_value() is None,
        )

        # Use the replacer.
        normalized = tuple(replacer.normalize(population))

    """

    def __init__(
        self,
        get_replacement,
        filter=lambda population, record: True,
    ):
        """
        Initialize a :class:`.ReplaceFitness` instance.

        Parameters
        ----------
        get_replacement : :class:`callable`
            Takes a single parameter, a :class:`tuple` of
            :class:`.MoleculeRecord` instances, which
            needs to be normalized. The entire population passed to
            :meth:`.normalize` is passed to this parameter, regardless
            of what is passed to the `filter` parameter. The
            :class:`callable` returns the value which is to be given
            to the normalized records.

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

        self._get_replacement = get_replacement
        self._filter = filter

    def normalize(self, population):
        replacement = self._get_replacement(population)
        for record in population:
            if self._filter(population, record):
                yield record.with_fitness_value(replacement)
            else:
                yield record
