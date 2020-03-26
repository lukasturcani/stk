"""
RAM Molecule Cache
==================

"""

from stk.molecular import InchiKey
from .molecule import MoleculeCache


class RamMoleculeCache(MoleculeCache):
    """
    Cache molecules from RAM.

    Most of the time you are better off using
    :class:`.MongDbMoleculeCache`. The :class:`.MongoDbMoleculeCache`
    includes a RAM based cache to prevent repeated reading and writing
    to disk, but the cache has a fixed size, so it will not grow
    indefinitely as more molecules are added, like
    :class:`.RamMoleculeCache` does. With
    :class:`.MongoDbMoleculeCache` if the
    RAM cache overflows its fixed size, the least recently used
    molecules are placed onto the hard disk, which is searched as a
    backup, if a molecule is not found in the RAM cache. Most of the
    time this will mean that :class:`.MongoDbMoleculeCache` has the
    same effective performance as :class:`.RamMoleculeCache`, but
    does not consume ever increasing RAM resources, and you end up
    with a proper molecular database when you are done using it,
    which can be used in future work. Note that for
    :class:`.MongoDbMoleculeCache` all molecules are written to the
    disk the first time they are put into the cache, so all molecules
    get into permanent storage, even if they remain in the RAM cache
    the entire time :class:`.MongoDbMoleculeCache` is in use, or if
    :class:`.MongoDbMoleculeCache` never overflows its fixed RAM size.

    """

    def __init__(self, key_maker=InchiKey()):
        """
        Initialize a :class:`.RamMoleculeCache` instance.

        Parameters
        ----------
        key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to create the key under which molecules get stored
            in the cache.

        """

        self._key_maker = key_maker
        self._cache = {}

    def put(self, molecule):
        self._cache[self._key_maker.get_key(molecule)] = molecule

    def get(self, key):
        return self._cache[key]
