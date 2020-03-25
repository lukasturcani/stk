"""
Constructed Molecule Cache
==========================

#. :class:`.MemoryConstructedMoleculeCache`
#. :class:`.MongoDbConstructedMoleculeCache`

"""


class ConstructedMoleculeCache:
    """
    Abstract base class for caching constructed molecules.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`constructed_molecule_cache \
    <~.caches.constructed_molecule.constructed_molecule>`, can
    serve as good examples.

    """

    def put(self, molecule):
        """
        Put `molecule` into the cache.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The molecule to place into the cache.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()

    def get(self, key, building_blocks, default=None):
        """
        Get the molecule with `key` from the cache.

        Parameters
        ----------
        key : :class:`object`
            The key of a molecule, which is to be returned from the
            cache.

        building_blocks : :class:`tuple` of :class:`.Molecule`
            The building blocks of the returned
            :class:`.ConstructedMolecule`. Their order should be equal
            their order in
            :meth:`.ConstructedMolecule.get_building_blocks`.

        default : :class:`.ConstructedMolecule`, optional
            If `key` is not found in the cache, return this molecule
            instead. If ``None``, an error will be raised if `key`
            is not found in the cache.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The molecule held in the cache under `key`.

        Raises
        ------
        :class:`KeyError`
            If `key` is not found in the cache and `default` is
            ``None``.

        """

        raise NotImplementedError()
