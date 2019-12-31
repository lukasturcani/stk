from .. import FunctionalGroup_


class Dibromo(FunctionalGroup_):
    """
    Represents a dibromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine1][atom1][atom2][bromine2]``.

    """

    def __init__(
        self,
        bromine1,
        atom1,
        bromine2,
        atom2,
        bonders,
        deleters,
    ):
        self._atom1 = atom1
        self._bromine1 = bromine1
        self._atom2 = atom2
        self._bromine2 = bromine2
        atoms = (bromine1, atom1, bromine2, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_bromine1(self):
        return self._bromine1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_bromine2(self):
        return self._bromine2.clone()
