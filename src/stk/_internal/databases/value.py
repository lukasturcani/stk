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

    *Iterating Through Entries in the Database*

    The :meth:`.MoleculeDatabase.get_all` and
    :meth:`.ConstructedMoleculeDatabase.get_all` methods  can be used
    to iterate through all of the database's entries, which can then
    be used to access values in the :class:`.ValueDatabase`

    .. testsetup:: iterating-through-entries-in-the-database

        import stk

        # Change the database used, so that when a developer
        # runs the doctests locally, their "stk" database is not
        # contaminated.
        _test_database = '_stk_doctest_database'

        _old_value_init = stk.ValueMongoDb
        stk.ValueMongoDb = lambda mongo_client, collection: (
            _old_value_init(
                mongo_client=mongo_client,
                database=_test_database,
                collection=collection,
            )
        )

        _old_molecule_init = stk.MoleculeMongoDb
        stk.MoleculeMongoDb = lambda mongo_client: _old_molecule_init(
            mongo_client=mongo_client,
            database=_test_database,
        )

        # Change the database MongoClient will connect to.

        import os
        import pymongo

        _mongo_client = pymongo.MongoClient
        _mongodb_uri = os.environ.get(
            'MONGODB_URI',
            'mongodb://localhost:27017/'
        )
        pymongo.MongoClient = lambda: _mongo_client(_mongodb_uri)

    .. testcode:: iterating-through-entries-in-the-database

        import pymongo

        client = pymongo.MongoClient()
        molecule_db = stk.MoleculeMongoDb(client)

        value_db = stk.ValueMongoDb(
            mongo_client=client,
            collection='atom_counts',
        )

        molecule = stk.BuildingBlock('BrBr')
        molecule_db.put(molecule)
        value_db.put(molecule, molecule.get_num_atoms())


        for molecule in molecule_db.get_all():
            try:
                value = value_db.get(molecule)
                print(value)
            except KeyError:
                # In case molecule is not in value_db.
                pass

    .. testoutput:: iterating-through-entries-in-the-database

        2

    .. testcleanup:: iterating-through-entries-in-the-database

        stk.ValueMongoDb = _old_value_init
        stk.MoleculeMongoDb = _old_molecule_init
        pymongo.MongoClient().drop_database(_test_database)
        pymongo.MongoClient = _mongo_client

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
