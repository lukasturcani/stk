from .fitness_normalizer import FitnessNormalizer


class Sum(FitnessNormalizer):
    """
    Sums the values in a :class:`list`.

    Examples
    --------
    .. code-block:: python

        # Create the normalizer.
        sum_normalizer = stk.Sum()

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = sum_normalizer.normalize(pop)

        # normalized is
        # {mol1: 6, mol2: 15, mol3: 24}

    """

    def __init__(self, filter=lambda population, mol: True):
        """
        Initialize a :class:`.Sum` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a two parameters, first is a :class:`.EAPopulation`
            and the second is a :class:`.Molecule`, and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.
            The :class:`.EAPopulation` on which :meth:`normalize` is
            called is passed as the first argument while the second
            argument will be passed every :class:`.Molecule` in it.

        """

        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        for mol in filtered:
            yield mol, sum(fitness_values[mol])
