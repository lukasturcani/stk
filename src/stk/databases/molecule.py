"""
Molecule Database
=================

.. toctree::
    :maxdepth: 2

    Molecule MongoDB <stk.databases.mongo_db.molecule>

"""


class MoleculeDatabase:
    """
    An abstract base class for storing and retrieving molecules.

    See Also
    --------
    :class:`.ConstructedMoleculeDatabase`
        If you need to store and retrieve
        :class:`.ConstructedMolecule` instances. You can put \
        :class:`.ConstructedMolecule` instances into a \
        :class:`.MoleculeDatabase`, however, you will only be able to \
        retrieve them as plain :class:`.Molecule` instances.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`molecule_database <~.databases.molecule>`, can serve as
    good examples.

    """

    def put(self, molecule):
        """
        Put `molecule` into the database.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to place into the database.

        Returns
        -------
        None : :class:`NoneType`

        """

        raise NotImplementedError()

    def get(self, key):
        """
        Get the molecule with `key` from the database.

        Parameters
        ----------
        key : :class:`object`
            The key of a molecule, which is to be returned from the
            database.

        Returns
        -------
        :class:`.Molecule`
            The molecule held in the database under `key`.

        Raises
        ------
        :class:`KeyError`
            If `key` is not found in the database.

        """

        raise NotImplementedError()

    def get_all(self):
        """
        Get all molecules in the database.

        Yields
        ------
        :class:`.Molecule`
            A molecule in the database.

        """

        raise NotImplementedError()
