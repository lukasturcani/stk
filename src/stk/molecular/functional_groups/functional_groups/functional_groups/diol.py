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
        atoms = (atom1, oxygen1, hydrogen1, atom2, oxygen2, hydrogen2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1

    def get_oxygen1(self):
        return self._oxygen1

    def get_hydrogen1(self):
        return self._hydrogen1

    def get_atom2(self):
        return self._atom2

    def get_oxygen2(self):
        return self._oxygen2

    def get_hydrogen2(self):
        return self._hydrogen2

    def clone(self):
        clone = super().clone()
        clone._atom1 = self._atom1
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._atom2 = self._atom2
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._atom1}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._atom2}, {self._oxygen2}, {self._hydrogen2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
