from .. import FunctionalGroup_


class Dibromo(FunctionalGroup_):
    """
    Represents a dibromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine1][atom1][atom2][bromine2]``.

    """

    def __init__(
        self,
        bromine1,
        atom1,
        bromine2,
        atom2,
        bonders,
        deleters,
    ):
        self._bromine1 = bromine1
        self._atom1 = atom1
        self._bromine2 = bromine2
        self._atom2 = atom2
        atoms = (bromine1, atom1, bromine2, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1

    def get_bromine1(self):
        return self._bromine1

    def get_atom2(self):
        return self._atom2

    def get_bromine2(self):
        return self._bromine2

    def clone(self):
        clone = super().clone()
        clone._atom1 = self._atom1
        clone._bromine1 = self._bromine1
        clone._atom2 = self._atom2
        clone._bromine2 = self._bromine2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._bromine1 = atom_map.get(
            self._bromine1.get_id(),
            self._bromine1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._bromine2 = atom_map.get(
            self._bromine2.get_id(),
            self._bromine2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine1}, {self._atom1}, {self._bromine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
