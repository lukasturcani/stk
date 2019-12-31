from .. import FunctionalGroup_


class Thioacid(FunctionalGroup_):
    """
    Represents a thioacid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[sulfur][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        sulfur,
        hydrogen,
        atom,
        bonders,
        deleters
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_sulfur(self):
        return self._sulfur.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()
