"""
RAM Molecule Value Cache
========================

"""

from stk.molecular import InchiKey
from .molecule import MoleculeValueCache


class RamMoleculeValueCache(MoleculeValueCache):
    """
    Cache molecular property values in RAM.

    Most of the time you are better off using
    :class:`.MongoDbMoleculeValueCache`. The
    :class:`.MongoDbMoleculeValueCache`
    includes a RAM based cache to prevent repeated reading and writing
    to disk, but the cache has a fixed size, so it will not grow
    indefinitely as more values are added, like
    :class:`.RamMoleculeValueCache` does. With
    :class:`.MongoDbMoleculeValueCache` if the
    RAM cache overflows its fixed size, the least recently used
    values are placed onto the hard disk, which is searched as a
    backup, if a value is not found in the RAM cache. Most of the
    time this will mean that :class:`.MongoDbMoleculeValueCache` has
    the same effective performance as :class:`.RamMoleculeValueCache`,
    but does not consume ever increasing RAM resources, and you end up
    with a proper database when you are done using it,
    which can be used in future work. Note that for
    :class:`.MongoDbMoleculeValueCache` all values are written to the
    disk the first time they are put into the cache, so all values
    get into permanent storage, even if they remain in the RAM cache
    the entire time :class:`.MongoDbMoleculeValueCache` is in use, or
    if :class:`.MongoDbMoleculeValueCache` never overflows its fixed
    RAM size.

    See Also
    --------
    :class:`.RamConstructedMoleculeValueCache`
        If you only need to use :class:`.MoleculeKeyMaker`, use a
        :class:`.RamMoleculeValueCache` instead of \
        :class:`.RamConstructedMoleculeValueCache`, \
        even if you are storing values associated with a \
        :class:`.ConstructedMolecule`.

    Examples
    --------
    You want to cache property values

    .. code-block:: python

        import stk

        molecule = stk.BuildingBlock('BrCCBr')
        ram_cache = stk.RamMoleculeValueCache()
        # Place a value into the cache.
        ram_cache.put(molecule, molecule.get_num_atoms())
        # Retrieve the value from the cache.
        num_atoms = ram_cache.get(molecule)

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
        ram_cache.put(polymer, polymer.get_num_atoms())
        num_polymer_atoms = ram_cache.get(polymer)

    """

    def __init__(self, key_maker=InchiKey()):
        """
        Initialize a :class:`.RamMoleculeValueCache` instance.

        Parameters
        ----------
        key_maker : :class:`.MoleculeKeyMaker`
            Used to get the key, which is associated with the values
            placed in the cache. If two molecules have
            the same key, they will return the same value from
            the cache.

        """

        self._key_maker = key_maker
        self._cache = {}

    def get(self, molecule):
        return self._cache[self._key_maker.get_key(molecule)]

    def put(self, molecule, value):
        self._cache[self._key_maker.get_key(molecule)] = value