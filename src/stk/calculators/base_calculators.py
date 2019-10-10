"""
Base Calculators
================

#. :class:`.Calculator`
#. :class:`.MoleculeCalculator`
#. :class:`.EAOperation`

This module provides calculator classes which serve as abstract base
classes for other types for calculators. Note that calculators defined
here do not act as direct base classes for calculators which are used
by users. The ususal inheritance scheme is :class:`.Calculator` is
subclassed by :class:`CalculatorType` which is subclassed by
:class:`UserCalculatorType`. For example, :class:`.Calculator` is
subclassed by :class:`.Optimizer` and :class:`.Optimizer` is
subclassed by :class:`.ETKDG`. Only :class:`.ETKDG` is instantiated
and used by the user. Note that :class:`.Optimizer` does not
subclass :class:`.Calculator` directly, it is subclassed via
:class:`.MolecularCalculator`, this is an implementation detail of the
:class:`.Optimizer` class.

:class:`.Calculator` simply serves as a common abstract base class
for every other ``stk`` calculator. The direct subclass of
:class:`.Calculator` is an abstract base class which defines a new
type of calculator for ``stk``, for example :class:`.Optimizer`
or :class:`.EnergyCalculator`. These define an interface for the new
type of calculation. Finally the subclasses of :class:`.Optimizer` or
:class:`EnergyCalculator` implement the calculation. For example
:class:`.ETKDG` or :class:`.MMFF``or different implementations of an
optimization, which the user can use.

For some of the abstract base classes provided here, an implementation
is also provided. If the implementation is inherited, the abstract base
class it implements should also be inherited, explicitly and directly.
This is because inheriting the implementation class is purely an
implementation detail that should be treated as invisible to the user
of the class.

"""


class Calculator:
    """
    Abstract base class for all calculators.

    """


class MoleculeCalculator(Calculator):
    """
    Abstract base class for calculators used on single molecules.

    """

    def _cache_result(self, fn, mol):
        """
        Apply `fn` to `mol` and cache the result.

        Parameters
        ----------
        fn : :class:`callable`
            Takes a :class:`.Molecule` as a parameter and returns
            some result, which is cached, depending on if
            caching is turned on.

        mol : :class:`.Molecule`
            The molecule onto which `fn` is applied.

        Returns
        -------
        :class:`object`
            The value returned by ``fn(mol)``.

        """

        if self.is_caching() and self.is_in_cache(mol):
            return self.get_cached_value(mol)

        value = fn(mol)

        if self.is_caching():
            self.add_to_cache(mol, value)

        return value

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

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def is_caching(self):
        """
        ``True`` if the calculator has caching turned on.

        Returns
        -------
        :class:`bool`
            ``True`` if the calculator has caching turned on.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

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

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

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

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def get_cached_value(self, mol):
        """
        Return the value stored in the cache for `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose cached value is to be returned.

        Returns
        -------
        :class:`object`
            The cached value.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class _MoleculeCalculator(MoleculeCalculator):
    """
    Implements the :class:`.MoleculeCalculator` interface.

    """

    def __init__(self, use_cache=False):
        """
        Initialize the :class:`.Calculator`.

        Parameters
        ----------
        use_cache : :class:`bool`, optional
            If ``True``, a calculation will not be performed on the
            same molecule twice.

        """

        self._cache = {}
        self._use_cache = use_cache

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

    def get_cached_value(self, mol):
        """
        Return the value stored in the cache for `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose cached value is to be returned.

        Returns
        -------
        :class:`object`
            The cached value.

        """

        return self._cache[mol]
