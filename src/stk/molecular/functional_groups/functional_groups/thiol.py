from .. import FunctionalGroup_


class Thiol(FunctionalGroup_):
    """
    Represents a thiol functional group.


    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][sulfur][hydrogen]``.

    """

    def __init__(self, sulfur, hydrogen, atom, bonders, deleters):
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_sulfur(self):
        return self._sulfur.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()
