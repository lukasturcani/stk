from .. import FunctionalGroup_


class Diol(FunctionalGroup_):
    """
    Represents a diol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[hydrogen1][oxygen1][atom1][atom2][oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        atom1,
        oxygen1,
        hydrogen1,
        atom2,
        oxygen2,
        hydrogen2,
        bonders,
        deleters,
    ):
        self._atom1 = atom1
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._atom2 = atom2
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        atoms = (
            atom1,
            oxygen1,
            hydrogen1,
            atom2,
            oxygen2,
            hydrogen2,
        )
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._atom1,
            self._oxygen1,
            self._hydrogen1,
            self._atom2,
            self._oxygen2,
            self._hydrogen2,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._atom1 = atom_map[self._atom1.id]
        clone._oxygen1 = atom_map[self._oxygen1.id]
        clone._hydrogen1 = atom_map[self._hydrogen1.id]
        clone._atom2 = atom_map[self._atom2.id]
        clone._oxygen2 = atom_map[self._oxygen2.id]
        clone._hydrogen2 = atom_map[self._hydrogen2.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._atom1}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._atom2}, {self._oxygen2}, {self._hydrogen2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
