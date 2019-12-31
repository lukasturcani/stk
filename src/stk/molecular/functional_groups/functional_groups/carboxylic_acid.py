from .. import FunctionalGroup_


class CarboxylicAcid(FunctionalGroup_):
    """
    Represents a carboxylic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen1])[oxygen2][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen1,
        oxygen2,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen1 = oxygen1
        self._oxygen2 = oxygen2
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen1, oxygen2, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()
