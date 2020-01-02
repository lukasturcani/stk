from .. import FunctionalGroup_


class BoronicAcid(FunctionalGroup_):
    """
    Represents a boronic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][boron]([oxygen1][hydrogen1])[oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        boron,
        oxygen1,
        hydrogen1,
        oxygen2,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._boron = boron
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_boron(self):
        return self._boron

    def get_oxygen1(self):
        return self._oxygen1

    def get_hydrogen1(self):
        return self._hydrogen1

    def get_oxygen2(self):
        return self._oxygen2

    def get_hydrogen2(self):
        return self._hydrogen2

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._boron = self._boron
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._boron = atom_map.get(
            self._boron.get_id(),
            self._boron,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._boron}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._oxygen2}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
