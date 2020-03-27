class CaseData:
    """
    A test case.

    Attributes
    ----------
    cache : class:`.MoleculeValueCache` or \
            :class:`.ConstructedMoleculeValueCache`
        The cache to test.

    molecule : :class:`.Molecule`
        The molecule to test.

    value : :class:`object`
        The value to put into the cache.

    """

    def __init__(self, cache, molecule, value):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        cache : class:`.MoleculeValueCache` or \
                :class:`.ConstructedMoleculeValueCache`
            The cache to test.

        molecule : :class:`.Molecule`
            The molecule to test.

        value : :class:`object`
            The value to put into the cache.

        """

        self.cache = cache
        self.molecule = molecule
        self.value = value
