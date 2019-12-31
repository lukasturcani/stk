from .. import FunctionalGroup_


class Difluoro(FunctionalGroup_):
    """
    Represents a difluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine1][atom1][atom2][fluorine2)``.

    """

    def __init__(
        self,
        atom1,
        fluorine1,
        atom2,
        fluorine2,
        bonders,
        deleters,
    ):
        self._atom1 = atom1
        self._fluorine1 = fluorine1
        self._atom2 = atom2
        self._fluorine2 = fluorine2
        atoms = (atom1, fluorine1, atom2, fluorine2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_fluorine1(self):
        return self._fluorine1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_fluorine2(self):
        return self._fluorine2.clone()
