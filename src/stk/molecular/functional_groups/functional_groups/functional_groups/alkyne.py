from .. import FunctionalGroup_


class Alkyne(FunctionalGroup_):
    """
    Represents an alkyne functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon1]([atom1])#[carbon2][atom2]``.

    """

    def __init__(
        self,
        carbon1,
        atom1,
        carbon2,
        atom2,
        bonders,
        deleters,
    ):
        self._carbon1 = carbon1
        self._atom1 = atom1
        self._carbon2 = carbon2
        self._atom2 = atom2
        atoms = (carbon1, atom1, carbon2, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1

    def get_carbon1(self):
        return self._carbon1

    def get_carbon2(self):
        return self._carbon2

    def get_atom2(self):
        return self._atom2

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._atom1 = self._atom1
        clone._carbon2 = self._carbon2
        clone._atom2 = self._atom2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._atom1}, {self._carbon2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
