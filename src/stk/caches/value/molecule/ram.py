"""
RAM Molecule Value Cache
========================

"""

from .molecule import MoleculeValueCache


class RamMoleculeValueCache(MoleculeValueCache):
    """
    Cache molecular property values in RAM.

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

