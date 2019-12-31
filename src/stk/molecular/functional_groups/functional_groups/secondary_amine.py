from .. import FunctionalGroup_


class SecondaryAmine(FunctionalGroup_):
    """
    Represents a secondary amine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom1][nitrogen]([hydrogen])[atom2]``.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen,
        atom1,
        atom2,
        bonders,
        deleters,
    ):
        self._nitrogen = nitrogen
        self._hydrogen = hydrogen
        self._atom1 = atom1
        self._atom2 = atom2
        atoms = (nitrogen, hydrogen, atom1, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom1(self):
        return self._atom1.clone()

    def get_atom2(self):
        return self._atom2.clone()
