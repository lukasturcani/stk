from .. import FunctionalGroup_


class Diol(FunctionalGroup_):
    """
    Represents a diol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[hydrogen1][oxygen1][atom1][atom2][oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        atom1,
        hydrogen1,
        oxygen1,
        atom2,
        oxygen2,
        hydrogen2,
        bonders,
        deleters,
    ):
        self._atom1 = atom1
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._atom2 = atom2
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        atoms = (
            atom1,
            oxygen1,
            hydrogen1,
            atom2,
            oxygen2,
            hydrogen2,
        )
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()
