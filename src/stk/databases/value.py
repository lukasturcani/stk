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

    *Iterating Through Entries in Database*

    The :meth:`.get_all` method of a molecule database instance can
    be used to iterate through entries, the keys of which can be used
    to access values.

    .. code-block:: python

        client = pymongo.MongoClient()
        molecule_db = stk.MoleculeMongoDb(client)

        value_db = stk.ValueMongoDb(
            mongo_client=client,
            collection='atom_counts',
        )

        for molecule in molecule_db.get_all():
            try:
                value = value_db.get(molecule)
            except KeyError:
                # In case molecule is not in value_db.
                pass

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
