from .. import FunctionalGroup_


class Alcohol(FunctionalGroup_):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(self, oxygen, hydrogen, atom, bonders, deleters):
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (oxygen, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()
