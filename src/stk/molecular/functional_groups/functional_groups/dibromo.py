from .. import FunctionalGroup_


class Dibromo(FunctionalGroup_):
    """
    Represents a dibromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom3][atom1]([bromine1])[atom2]([bromine2])[atom4]``.

    """

    def __init__(
        self,
        atom1,
        bromine1,
        atom2,
        bromine2,
        atom3,
        atom4,
        bonders,
        deleters,
    ):
        self._atom1 = atom1
        self._bromine1 = bromine1
        self._atom2 = atom2
        self._bromine2 = bromine2
        self._atom3 = atom3
        self._atom4 = atom4
        atoms = (atom1, bromine1, atom2, bromine2, atom3, atom4)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_bromine1(self):
        return self._bromine1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_bromine2(self):
        return self._bromine2.clone()

    def get_atom3(self):
        return self._atom3.clone()

    def get_atom4(self):
        return self._atom4.clone()
