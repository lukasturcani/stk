from .. import FunctionalGroup_


class Amide(FunctionalGroup_):
    """
    Represents an amide functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        nitrogen,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._carbon,
            self._oxygen,
            self._nitrogen,
            self._hydrogen1,
            self._hydrogen2,
            self._atom,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.id]
        clone._oxygen = atom_map[self._oxygen.id]
        clone._nitrogen = atom_map[self._nitrogen.id]
        clone._hydrogen1 = atom_map[self._hydrogen1.id]
        clone._hydrogen2 = atom_map[self._hydrogen2.id]
        clone._atom = atom_map[self._atom.id]
        return clone
