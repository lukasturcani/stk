"""
Value MongoDB
=============

"""

from functools import lru_cache

from stk.molecular import InchiKey

from ..value import ValueDatabase
from .utilities import HashableDict


class ValueMongoDb(ValueDatabase):
    """
    Use MongoDB to store and retrieve molecular property values.

    Examples
    --------
    See also examples in :class:`.ValueDatabase`.

    *Storing Molecular Properties in a Database*

    You want to store property values in a database.

    .. testsetup:: storing-molecular-properties-in-a-database

        import stk

        # Change the database used, so that when a developer
        # runs the doctests locally, their "stk" database is not
        # contaminated.
        _test_database = '_stk_doctest_database'
        _old_init = stk.ValueMongoDb
        stk.ValueMongoDb = lambda mongo_client, collection: (
            _old_init(
                mongo_client=mongo_client,
                database=_test_database,
                collection=collection,
            )
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

    .. testcode:: storing-molecular-properties-in-a-database

        import stk
        import pymongo

        # Connect to a MongoDB. This example connects to a local
        # MongoDB, but you can connect to a remote DB too with
        # MongoClient() - read the documentation for pymongo to see how
        # to do that.
        client = pymongo.MongoClient()
        db = stk.ValueMongoDb(
            mongo_client=client,
            collection='atom_counts',
        )

        molecule = stk.BuildingBlock('BrCCBr')
        # Add the value to the database.
        db.put(molecule, molecule.get_num_atoms())
        # Retrieve the value from the database.
        num_atoms = db.get(molecule)

        # Works with constructed molecules too.
        polymer = stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
                ),
                repeating_unit='A',
                num_repeating_units=2,
            ),
        )
        db.put(polymer, polymer.get_num_atoms())
        num_polymer_atoms = db.get(polymer)

    .. testcode:: storing-molecular-properties-in-a-database
        :hide:

        assert num_polymer_atoms == polymer.get_num_atoms()

    .. testcleanup:: storing-molecular-properties-in-a-database

        stk.ValueMongoDb = _old_init
        pymongo.MongoClient().drop_database(_test_database)
        pymongo.MongoClient = _mongo_client

    """

    def __init__(
        self,
        mongo_client,
        collection,
        database='stk',
        key_makers=(InchiKey(), ),
        put_lru_cache_size=128,
        get_lru_cache_size=128,
        indices=('InChIKey', ),
    ):
        """
        Initialize a :class:`.ValueMongoDb` instance.

        Parameters
        ----------
        mongo_client : :class:`pymongo.MongoClient`
            The database client.

        collection : :class:`str`
            The name of the MongoDB collection used for storing the
            property values.

        database : :class:`str`, optional
            The name of the MongoDB database used for storing the
            property values.

        key_makers : :class:`tuple` of :class:`.MoleculeKeyMaker`
            Used to make the keys of molecules, which the values
            are associated with. If two molecules have the same
            key, they will return the same value from the database.

        put_lru_cache_size : :class:`int`, optional
            A RAM-based least recently used cache is used to avoid
            writing to the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        get_lru_cache_size : :class:`int`, optional
            A RAM-based least recently used cache is used to avoid
            reading from the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        indices : :class:`tuple` of :class:`str`, optional
            The names of molecule keys, on which an index should be
            created, in order to minimize lookup time.

        """

        self._values = mongo_client[database][collection]
        self._key_makers = key_makers
        self._put = lru_cache(maxsize=put_lru_cache_size)(self._put)
        self._get = lru_cache(maxsize=get_lru_cache_size)(self._get)

        index_information = self._values.index_information()
        if 'v_1' not in index_information:
            self._values.create_index('v')

        for index in indices:
            # Do not create the same index twice.
            if f'{index}_1' not in index_information:
                self._values.create_index(index)

    def put(self, molecule, value):
        json = {'v': value}
        for key_maker in self._key_makers:
            json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )
        # lru_cache requires that the parameters to the cached function
        # are hashable objects.
        return self._put(HashableDict(json))

    def _put(self, json):
        keys = dict(json)
        keys.pop('v')

        query = {'$or': []}
        for key, value in keys.items():
            query['$or'].append({key: value})

        self._values.update_many(
            filter=query,
            update={
                '$set': json
            },
            upsert=True,
        )

    def get(self, molecule):

        def make_dict(key_maker):
            return HashableDict({
                key_maker.get_key_name():
                key_maker.get_key(molecule)
            })

        key = {'$or': tuple(map(make_dict, self._key_makers))}
        # lru_cache requires that the parameters to the cached function
        # are hashable objects.
        return self._get(HashableDict(key))

    def _get(self, key):
        value = self._values.find_one(key)
        if value is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        return value['v']
