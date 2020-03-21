"""
Molecule Key
============

"""


class MoleculeKey:
    """
    An abstract base class for molecular keys.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is purely for convenience when
    implementing subclasses. The implemented public methods are
    simply default implementations, which can be safely ignored or
    overridden, when implementing subclasses.



    """

    def __init__(self, name, key):
        self._name = name
        self._key = key

    def get_name(self):
        return self._name

    def get_key(self, molecule):
        return self._key(molecule)
