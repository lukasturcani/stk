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
        atom_map = {
            bromine1.id: bromine1.clone(),
            atom1.id: atom1.clone(),
            bromine2.id: bromine2.clone(),
            atom2.id: atom2.clone(),
        }
        self._atom1 = atom_map[atom1.id]
        self._bromine1 = atom_map[bromine1.id]
        self._atom2 = atom_map[atom2.id]
        self._bromine2 = atom_map[bromine2.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.id] for a in bonders),
            deleters=tuple(atom_map[a.id] for a in deleters),
        )

    def get_atom1(self):
        return self._atom1.clone()

    def get_bromine1(self):
        return self._bromine1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def get_bromine2(self):
        return self._bromine2.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._atom1, self._bromine1, self._atom2, self._bromine2
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._atom1 = atom_map[self._atom1.id]
        clone._bromine1 = atom_map[self._bromine1.id]
        clone._atom2 = atom_map[self._atom2.id]
        clone._bromine2 = atom_map[self._bromine2.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine1}, {self._atom1}, {self._bromine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
