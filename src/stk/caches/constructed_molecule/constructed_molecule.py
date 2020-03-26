"""
Constructed Molecule Cache
==========================

#. :class:`.RamConstructedMoleculeCache`
#. :class:`.MongoDbConstructedMoleculeCache`

"""

from ..molecule import MoleculeCache


class ConstructedMoleculeCache(MoleculeCache):
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
        :class:`.ConstructedMolecule`
            The molecule held in the cache under `key`.

        Raises
        ------
        :class:`KeyError`
            If `key` is not found in the cache.

        """

        raise NotImplementedError()
