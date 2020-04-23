class CaseData:
    """
    A test case.

    Attributes
    ----------
    fitness_normalizer : :class:`.FitnessNormalizer`
        The fitness normalizer to test.

    population : :class:`tuple` of :class:`.MoleculeRecord`
        The population which is normalized.

    normalized : :class:`tuple` of :class:`.MoleculeRecord`
        The normalized :attr:`.population`.

    """

    def __init__(self, fitness_normalizer, population, normalized):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        fitness_normalizer : :class:`.FitnessNormalizer`
            The fitness normalizer to test.

        population : :class:`tuple` of :class:`.MoleculeRecord`
            The population which is normalized.

        normalized : :class:`tuple` of :class:`.MoleculeRecord`
            The normalized `population`.

        """

        self.fitness_normalizer = fitness_normalizer
        self.population = population
        self.normalized = normalized
