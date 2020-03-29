from .fitness_normalizer import FitnessNormalizer


class Power(FitnessNormalizer):
    """
    Raises fitness values to some power.

    This works for cases where the fitness value is single
    :class:`float` and where it is :class:`list` of :class:`float`.

    Examples
    --------
    Raising a fitness value by some power

    .. code-block:: python

        import stk

        pop = stk.Population(...)
        # Assume this returns {mol1: 1, mol2: 2, mol3: 3}.
        pop.get_fitness_values()

        # Create the normalizer.
        power = stk.Power(2)

        # Normalize the fitness values.
        normalized = power.normalize(pop)

        # normalized is {mol1: 1, mol2: 4, mol3: 9}.


    Raising vector valued fitness values by some power

    .. code-block:: python

        # Create the normalizer.
        power = stk.Power(2)

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = power.normalize(pop)

        # normalized is
        # {mol1: [1, 4, 9], mol2: [16, 25, 36], mol3: [49, 64, 81]}.


    Raising vector valued fitness values by different powers

    .. code-block:: python

        # Create the normalizer.
        power = stk.Power([1, 2, 3])

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.

        # Normalize the fitness values.
        normalized = power.normalize(pop)

        # normalized is
        # {mol1: [1, 4, 27], mol2: [4, 25, 216], mol3: [7, 64, 729]}.

    """

    def __init__(self, power, filter=lambda population, mol: True):
        """
        Initialize a :class:`Power` instance.

        Parameters
        ----------
        power : :class:`float` or :class:`list` of :class:`float`
            The power to raise each :attr:`fitness` value to. Can be
            a single number or multiple numbers.

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

        self._power = power
        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        for mol in filtered:
            yield mol, np.float_power(fitness_values[mol], self._power)


