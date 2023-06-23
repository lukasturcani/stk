class FitnessNormalizer:
    """
    Abstract base class for fitness normalizers.

    A fitness normalizer takes an :class:`tuple` of
    :class:`.MoleculeRecord` instances and yields new
    :class:`.MoleculeRecord` instances, with normalized fitness values.
    The primary benefit of a normalizer vs a
    :class:`.FitnessCalculator` is that a :class:`.FitnessNormalizer`
    has access to all members in the population when it is calculating
    the normalized fitness value, whereas a :class:`.FitnessCalculator`
    does not.

    """

    def normalize(self, population):
        """
        Normalize the fitness values in `population`.

        Parameters
        ----------
        population : :class:`tuple` of :class:`.MoleculeRecord`
            The molecules which need to have their fitness values
            normalized.

        Yields
        ------
        :class:`.MoleculeRecord`
            A record with a normalized fitness value.

        """

        raise NotImplementedError()
