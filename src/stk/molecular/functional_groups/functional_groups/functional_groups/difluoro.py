from .. import FunctionalGroup_


class Difluoro(FunctionalGroup_):
    """
    Represents a difluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine1][atom1][atom2][fluorine2]``.

    """

    def __init__(
        self,
        fluorine1,
        atom1,
        fluorine2,
        atom2,
        bonders,
        deleters,
    ):
        self._fluorine1 = fluorine1
        self._atom1 = atom1
        self._fluorine2 = fluorine2
        self._atom2 = atom2
        atoms = (fluorine1, atom1, fluorine2, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_atom1(self):
        return self._atom1

    def get_fluorine1(self):
        return self._fluorine1

    def get_atom2(self):
        return self._atom2

    def get_fluorine2(self):
        return self._fluorine2

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._fluorine1 = atom_map.get(
            self._fluorine1.get_id(),
            self._fluorine1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._fluorine2 = atom_map.get(
            self._fluorine2.get_id(),
            self._fluorine2,
        )
        return clone

    def clone(self):
        clone = super().clone()
        clone._atom1 = self._atom1
        clone._fluorine1 = self._fluorine1
        clone._atom2 = self._atom2
        clone._fluorine2 = self._fluorine2
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine1}, {self._atom1}, {self._fluorine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
