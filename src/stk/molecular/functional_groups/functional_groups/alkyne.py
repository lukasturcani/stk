from .. import FunctionalGroup_


class Alkyne(FunctionalGroup_):
    """
    Represents an alkyne functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon1]([atom1])#[carbon2][atom2]``.

    """

    def __init__(
        self,
        carbon1,
        atom1,
        carbon2,
        atom2,
        bonders,
        deleters,
    ):
        self._carbon1 = carbon1
        self._atom1 = atom1
        self._carbon2 = carbon2
        self._atom2 = atom2
        atoms = (carbon1, atom1, carbon2, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_carbon1(self):
        return self._carbon1.clone()

    def get_carbon2(self):
        return self._carbon2.clone()

    def get_atom2(self):
        return self._atom2.clone()
