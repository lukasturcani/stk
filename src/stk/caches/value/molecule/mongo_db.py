"""
MongoDB Value Cache
===================

"""

from functools import lru_cache

from stk.molecular import InchiKey
from .value_cache import ValueCache
from ...utilities import HashableDict


class MongoDbMoleculeValueCache(ValueCache):
    """
    Use MongoDB to store and retrieve molecular property values.

    See Also
    --------
    :class:`.MongoDbConstructedMoleculeValueCache`
        If you need to use :class:`.ConstructedMoleculeKeyMaker` \
        instead of, or in addition to, :class:`.MoleculeKeyMaker`. \
        If you only need to use :class:`.MoleculeKeyMaker` use this \
        class, even if you are storing values associated with a \
        :class:`.ConstructedMolecule`.

    Examples
    -------
    Before using this class, make sure you have :mod:`pymongo` and
    that its working properly. I recommend reading at least the
    introductory and installation
    documentation of :mod:`pymongo` before using this class. Those
    docs can be found here__.

    __ https://api.mongodb.com/python/current/

    You want

    """

    def __init__(
        self,
        mongo_client,
        collection,
        database='stk',
        key_makers=(InchiKey(), ),
        lru_cache_size=128,
    ):
        """
        Initialize a :class:`.MongDbValueCache` instance.

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

        key_makers : :class:`tuple` of :class:`MoleculeKeyMaker`
            Used to make the keys of molecules, which the values
            are associated with.

        lru_cache_size : :class:`int`, optional
            A RAM-based least recently used cache is used to avoid
            reading and writing to the database repeatedly. This sets
            the number of values which fit into the LRU cache. If
            ``None``, the cache size will be unlimited.

        """

        self._values = mongo_client[database][collection]
        self._key_makers = key_makers
        self._put = lru_cache(maxsize=lru_cache_size)(self._put)
        self._get = lru_cache(maxsize=lru_cache_size)(self._get)

    def put(self, molecule, value):
        json = {'v': value}
        for key_maker in self._key_makers:
            json[key_maker.get_key_name()] = (
                key_maker.get_key(molecule)
            )

        return self._put(json)

    def _put(self, json):
        self._values.insert_one(json)

    def get(self, molecule):
        key = {
            '$or': [
                {key_maker.get_key_name(): key_maker.get_key(molecule)}
                for key_maker in self._key_makers
            ],
        }
        return self._get(HashableDict(key))

    def _get(self, key):
        value = self._values.find_one(key)
        if value is None:
            raise KeyError(
                'No molecule found in the database with a key of: '
                f'{key}'
            )
        return value
