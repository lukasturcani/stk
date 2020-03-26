class CaseData:
    """
    A test case.

    Attributes
    ----------
    cache : :class:`.MoleculeCache
        The cache to test.

    molecule : :class:`.Molecule`
        The molecule to put and get from the :attr:`.cache`.

    key : :class:`object`
        The key used to retrieve :attr:`.molecule` from the cache.

    """

    def __init__(self, cache, molecule, key):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        cache : :class:`.MoleculeCache
            The cache to test.

        molecule : :class:`.Molecule`
            The molecule to put and get from the `cache`.

        key : :class:`object`
            The key used to retrieve `molecule` from the cache.

        """

        self.cache = cache
        self.molecule = molecule
        self.key = key
