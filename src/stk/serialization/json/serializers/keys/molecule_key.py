class MoleculeKey:
    """

    """

    def __init__(self, name, key):
        self._name = name
        self._key = key

    def get_name(self):
        return self._name

    def get_key(self, molecule):
        return self._key(molecule)
