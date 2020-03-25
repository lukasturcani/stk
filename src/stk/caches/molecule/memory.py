from . molecule import MoleculeCache


class MemoryMoleculeCache(MoleculeCache):
    """

    """

    def __init__(self, key_maker):
        """

        """

        self._key_maker = key_maker
        self._cache = {}

    def put(self, molecule):
        self._cache[self._key_maker.get_key(molecule)] = molecule

    def get(self, key, default=None):
        if default is None:
            return self._cache[key]
        return self._cache.get(key, default)
