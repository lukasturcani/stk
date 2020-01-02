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
            atom1.get_id(): atom1.clone(),
            oxygen1.get_id(): oxygen1.clone(),
            hydrogen1.get_id(): hydrogen1.clone(),
            atom2.get_id(): atom2.clone(),
            oxygen2.get_id(): oxygen2.clone(),
            hydrogen2.get_id(): hydrogen2.clone(),
        }
        self._atom1 = atom_map[atom1.get_id()]
        self._oxygen1 = atom_map[oxygen1.get_id()]
        self._hydrogen1 = atom_map[hydrogen1.get_id()]
        self._atom2 = atom_map[atom2.get_id()]
        self._oxygen2 = atom_map[oxygen2.get_id()]
        self._hydrogen2 = atom_map[hydrogen2.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
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
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._atom1 = atom_map[self._atom1.get_id()]
        clone._oxygen1 = atom_map[self._oxygen1.get_id()]
        clone._hydrogen1 = atom_map[self._hydrogen1.get_id()]
        clone._atom2 = atom_map[self._atom2.get_id()]
        clone._oxygen2 = atom_map[self._oxygen2.get_id()]
        clone._hydrogen2 = atom_map[self._hydrogen2.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._atom1}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._atom2}, {self._oxygen2}, {self._hydrogen2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
