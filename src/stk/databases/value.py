"""
Value Database
==============

.. toctree::
    :maxdepth: 2

    Value MongoDB <stk.databases.mongo_db.value>

"""


class ValueDatabase:
    """
    Abstract base class for storing of molecular property values.

    Examples
    --------
    *Subclass Implementation*

    The source of any of the subclasses, listed in
    :mod:`molecule_value_database <.databases.value>`,
    can serve as good examples.

    """

    def put(self, molecule, value):
        """
        Put a value into the database.

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
        Get the stored value for `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule whose value is to be retrieved from the
            database.

        Returns
        -------
        :class:`object`
            The value associated with `molecule`.

        Raises
        ------
        :class:`KeyError`
            If `molecule` is not found in the database.

        """

        raise NotImplementedError()
