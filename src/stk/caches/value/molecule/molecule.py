"""
Molecule Value Cache
====================

#. :class:`.MongoDbValueCache`
#. :class:`.RamValueCache`

"""

from ..constructed_molecule import ConstructedMoleculeValueCache


class MoleculeValueCache(ConstructedMoleculeValueCache):
    """
    Abstract base class for caching of molecular property values.

    Note that a :class:`.MoleculeValueCache` can be used anywhere a
    :class:`.ConstructedMoleculeValueCache` is required. However, the
    opposite is not true. If something requires a
    :class:`.MoleculeValueCache` you cannot use a
    :class:`.ConstructedMoleculeValueCache` in its place.

    """

    def put(self, molecule, value):
        """
        Put a value into the cache.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule which is associated with the `value`.

        value : :class:`object`
            Some value associated with `molecule`.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()

    def get(self, molecule):
        """
        Get the cached value for `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule whose value is to be retrieved from the cache.

        Returns
        -------
        :class:`object`
            The value associated with `molecule`.

        Raises
        ------
        :class:`KeyError`
            If `molecule` is not found in the cache.

        """

        raise NotImplementedError()
