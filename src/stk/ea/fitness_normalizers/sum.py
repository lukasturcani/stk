"""
Sum
===

"""

from .fitness_normalizer import FitnessNormalizer


class Sum(FitnessNormalizer):
    """
    Sums fitness values .

    Examples
    --------
    *Combining Fitness Value Into a Single Value*

    .. code-block:: python

        # Create the normalizer.
        sum_normalizer = stk.Sum()

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = sum_normalizer.normalize(pop)

        # normalized is
        # {mol1: 6, mol2: 15, mol3: 24}

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
