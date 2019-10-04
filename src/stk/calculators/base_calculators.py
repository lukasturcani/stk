"""
Base Calculators
================

#. :class:`.Calculator`
#. :class:`.MoleculeCalculator`

This module provide calculator classes which serve as bases classes
for other types for calculators. Note that calculators defined here do
not act as base classes for calculators which are used by users. The
ususal inheritance scheme is :class:`.Calculator` is subclassed by
:class:`CalculatorType` which is subclassed by
:class:`UserCalculatorType`. For example, :class:`.Calculator` is
subclassed by :class:`.Optimizer`
(technically via :class:`.MolecularCalculator` but that is an
implementation detail which doesn't matter) and :class:`.Optimizer` is
subclassed by :class:`.ETKDG`. Only :class:`ETKDG` is instantiated
and used by the user.

So :class:`.Calculator` simply serves as a common base class for
every other ``stk`` calculators. The direct subclass of
:class:`.Calculator` is an abstract base class which defines a new
type of calculator for ``stk``, for example :class:`.Optimizer`
or :class:`.EnergyCalculator`. These define the type of calculation.
Finally the subclasses of :class:`.Optimizer` or
:class:`EnergyCalculator` implement the calculation. For example
:class:`.ETKDG` or :class:`.MMFF``or different implementations of an
optimization, which the user can use.

Of course there can be intermediate classes between any of these
subclasses to allow code re-use, :class:`.MoleculeCalculator` is an
example of such a class.

"""


class Calculator:
    """
    A base class for all calculators.

    """

    def __init__(self, **kwargs):
        """
        Initialize a :class:`.Calculator`.

        """

        return


class MoleculeCalculator(Calculator):
    """
    Base class for calculators which operate on single molecules.

    """

    def __init__(self, use_cache=False, **kwargs):
        """
        Initialize a :class:`.MoleculeCalculator`

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True``, the calculator will not perform a calculation
            on the same molecule twice.

        """

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

        self._cache[mol] = value
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


class EAOperation(Calculator):
    """
    Base class for :class:`.Crosser` and :class:`.Mutator`.

    """

    def __init__(self, use_cache=False, **kwargs):
        self._use_cache = use_cache
        super().__init__(use_cache=use_cache, **kwargs)

    def set_cache_use(self, use_cache):
        """
        Set use of the molecular cache on or off.

        Parameters
        ----------
        use_cache : :class:`bool`
            ``True`` if the molecular cache is to be used.

        Returns
        -------
        :class:`.EAOperation`
            The calculator.

        """

        self._use_cache = use_cache
        return self
