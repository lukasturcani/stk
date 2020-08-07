"""
Constructed Molecule Database
=============================

.. toctree::
    :maxdepth: 2

    Constructed Molecule MongoDB <\
stk.databases.mongo_db.constructed_molecule\
>

"""


class ConstructedMoleculeDatabase:
    """
    Abstract base class for storing constructed molecules.

    See Also
    --------
    :class:`.MoleculeDatabase`
        If you need to store and retrieve :class:`.Molecule`
        instances.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in
    :mod:`constructed_molecule_database \
    <~.databases.constructed_molecule>`, can
    serve as good examples.

    """

    def put(self, molecule):
        """
        Put `molecule` into the database.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
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
        :class:`.ConstructedMolecule`
            The molecule held in the database under `key`.

        Raises
        ------
        :class:`KeyError`
            If `key` is not found in the database.

        """

        raise NotImplementedError()

    def get_all(self):
        """
        Get all entries in the database.

        Yields
        ------
        :class:`.ConstructedMolecule`
            A molecule in the database.

        """

        raise NotImplementedError()
