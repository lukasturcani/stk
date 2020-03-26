"""
Constructed Molecule Value Cache
================================

#. :class:`.MongoDbMoleculeValueCache`
#. :class:`.MongoDbConstructedMoleculeValueCache`
#. :class:`.MoleculeValueCache`
#. :class:`.RamConstructedMoleculeValueCache`
#. :class:`.RamMoleculeValueCache`

"""


class ConstructedMoleculeValueCache:
    """
    Abstract base class for caching of molecular property values.

    Note that a :class:`.MoleculeValueCache` can be used anywhere a
    :class:`.ConstructedMoleculeValueCache` is required. However, the
    opposite is not true. If something requires a
    :class:`.MoleculeValueCache` you cannot use a
    :class:`.ConstructedMoleculeValueCache` in its place.

    Examples
    --------
    *Subclass Implementation*

    The source of of any of the subclasses, listed in
    :mod:`constructed_molecule_value_cache \
    <.caches.value.constructed_molecule.constructed_molecule>`,
    can serve as good examples.

    """

    def put(self, molecule, value):
        """
        Put a value into the cache.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
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
        molecule : :class:`.ConstructedMolecule`
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
