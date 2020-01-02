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
        super().__init__(atoms, bonders, atoms)

    def get_carbon(self):
        return self._carbon

    def get_oxygen1(self):
        return self._oxygen1

    def get_oxygen2(self):
        return self._oxygen2

    def get_hydrogen(self):
        return self._hydrogen

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen1 = self._oxygen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen1}, {self._oxygen2}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
