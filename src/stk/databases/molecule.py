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

    # Keep empty __init__() to override ugly default docstring.
    def __init__(self):
        """"""
        return

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

    def put_many(self, molecules):
        """
        Put `molecules` into the database.

        Parameters
        ----------
        molecules : :class:`iterable` of :class:`.Molecule`
            The molecules to place into the database.

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
