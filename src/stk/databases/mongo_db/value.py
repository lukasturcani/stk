"""
Value MongoDB
=============

"""

from functools import lru_cache
import warnings

from stk.molecular import InchiKey
from .utilities import HashableDict
from ..value import ValueDatabase


class ValueMongoDb(ValueDatabase):
    """
    Use MongoDB to store and retrieve molecular property values.

    Examples
    --------
    *Storing Molecular Properties in a Database*

    You want to store property values in a database.

    .. code-block:: python

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
                num_repeating_units=2',
            ),
        )
        db.put(polymer, polymer.get_num_atoms())
        num_polymer_atoms = db.get(polymer)

    """

    def __init__(
        self,
        mongo_client,
        collection,
        database='stk',
        key_makers=(InchiKey(), ),
        lru_cache_size='',
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

        lru_cache_size : :class:`int`, optional
            This argument is deprecated and will be removed in any
            version of :mod:`stk` released on, or after, 01/01/21.
            Use the `put_lru_cache_size` and `get_lru_cache_size`
            arguments instead.

            A RAM-based least recently used cache is used to avoid
            reading and writing to the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

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

        if lru_cache_size != '':
            warnings.warn(
                'The lru_cache_size argument is deprecated and will '
                'be removed in any version of stk released on, or '
                'after, 01/01/21. Use the put_lru_cache_size and '
                'get_lru_cache_size arguments instead.',
                FutureWarning,
            )
            put_lru_cache_size = lru_cache_size
            get_lru_cache_size = lru_cache_size

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
