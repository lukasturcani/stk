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
        atom_map = {
            atom1.id: atom1.clone(),
            oxygen1.id: oxygen1.clone(),
            hydrogen1.id: hydrogen1.clone(),
            atom2.id: atom2.clone(),
            oxygen2.id: oxygen2.clone(),
            hydrogen2.id: hydrogen2.clone(),
        }
        self._atom1 = atom_map[atom1.id]
        self._oxygen1 = atom_map[oxygen1.id]
        self._hydrogen1 = atom_map[hydrogen1.id]
        self._atom2 = atom_map[atom2.id]
        self._oxygen2 = atom_map[oxygen2.id]
        self._hydrogen2 = atom_map[hydrogen2.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.id] for a in bonders),
            deleters=tuple(atom_map[a.id] for a in deleters),
        )

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
