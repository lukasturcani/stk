"""
Molecule Cache
===============

#. :class:`.RamMoleculeCache`
#. :class:`.MongoDbMoleculeCache`

"""


class MoleculeCache:
    """
    An abstract base class for caching molecules.

    See Also
    --------
    :class:`.ConstructedMoleculeCache`
        If you need to store and retrieve
        :class:`.ConstructedMolecule` instances. You can put \
        :class:`.ConstructedMolecule` instances into a \
        :class:`.MoleculeCache`, however, you will only be able to \
        retrieve them as plain :class:`.Molecule` instances.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`molecule_cache <~.caches.molecule.molecule>`, can serve as
    good examples.

    """

    def put(self, molecule):
        """
        Put `molecule` into the cache.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to place into the cache.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()

    def get(self, key):
        """
        Get the molecule with `key` from the cache.

        Parameters
        ----------
        key : :class:`object`
            The key of a molecule, which is to be returned from the
            cache.

        Returns
        -------
        :class:`.Molecule`
            The molecule held in the cache under `key`.

        Raises
        ------
        :class:`KeyError`
            If `key` is not found in the cache.

        """

        raise NotImplementedError()
