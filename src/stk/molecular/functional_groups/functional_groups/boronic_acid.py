from .. import FunctionalGroup_


class BoronicAcid(FunctionalGroup_):
    """
    Represents a boronic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][boron]([oxygen1][hydrogen1])[oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        boron,
        oxygen1,
        hydrogen1,
        oxygen2,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._boron = boron
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_boron(self):
        return self._boron.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()
