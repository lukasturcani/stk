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
        atom_map = {
            atom1.id: atom1.clone(),
            fluorine1.id: fluorine1.clone(),
            atom2.id: atom2.clone(),
            fluorine2.id: fluorine2.clone(),
        }
        self._atom1 = atom_map[atom1.id]
        self._fluorine1 = atom_map[fluorine1.id]
        self._atom2 = atom_map[atom2.id]
        self._fluorine2 = atom_map[fluorine2.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.id] for a in bonders),
            deleters=tuple(atom_map[a.id] for a in deleters),
        )

    def get_atom1(self):
        return self._atom1.clone()

    def get_fluorine1(self):
        return self._fluorine1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_fluorine2(self):
        return self._fluorine2.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._atom1, self._fluorine1, self._atom2, self._fluorine2
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._atom1 = atom_map[self._atom1.id]
        clone._fluorine1 = atom_map[self._fluorine1.id]
        clone._atom2 = atom_map[self._atom2.id]
        clone._fluorine2 = atom_map[self._fluorine2.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine1}, {self._atom1}, {self._fluorine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
