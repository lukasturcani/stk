"""
Base Calculators
================



"""


class Calculator:
    """
    A base class for all calculators.

    """

    def __init__(self, **kwargs):
        return


class MoleculeCalculator(Calculator):
    """
    Base class for calculators which operate on single molecules.

    """

    def __init__(self, use_cache=False, **kwargs):
        self._use_cache = use_cache
        self._cache = {}
        super().__init__(use_cache=use_cache, **kwargs)

    def set_cache_use(self, use_cache):
        """
        Set cache use on or off.

        Parameters
        ----------
        use_cache : :class:`bool`
            ``True`` if the cache is to be used.

        Returns
        -------
        :class:`.MoleculeCalculator`
            The calculator.

        """

        self._use_cache = use_cache
        return self

    def is_caching(self):
        """
        ``True`` if the calculator has caching turned on.

        Returns
        -------
        :class:`bool`
            ``True`` if the calculator has caching turned on.

        """

        return self._use_cache

    def add_to_cache(self, mol, value=None):
        """
        Add a molecule to the cache.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule to be added to the cache.

        value : class:`object`, optional
            The cached value associated with the molecule.

        Returns
        -------
        :class:`.MoleculeCalculator`
            The calculator.

        """

        self._cache.add(mol)
        return self

    def is_in_cache(self, mol):
        """
        Return ``True`` if `mol` is cached.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule being checked.

        Returns
        -------
        :class:`bool`
            ``True`` if `mol` is cached.

        """

        return mol in self._cache
