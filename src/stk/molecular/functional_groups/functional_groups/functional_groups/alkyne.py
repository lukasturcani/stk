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
        atom_map = {
            carbon1.id: carbon1.clone(),
            atom1.id: atom1.clone(),
            carbon2.id: carbon2.clone(),
            atom2.id: atom2.clone(),
        }
        self._carbon1 = atom_map[carbon1.id]
        self._atom1 = atom_map[atom1.id]
        self._carbon2 = atom_map[carbon2.id]
        self._atom2 = atom_map[atom2.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.id] for a in bonders),
            deleters=tuple(atom_map[a.id] for a in deleters),
        )

    def get_atom1(self):
        return self._atom1.clone()

    def get_carbon1(self):
        return self._carbon1.clone()

    def get_carbon2(self):
        return self._carbon2.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._carbon1, self._atom1, self._carbon2, self._atom2
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon1 = atom_map[self._carbon1.id]
        clone._atom1 = atom_map[self._atom1.id]
        clone._carbon2 = atom_map[self._carbon2.id]
        clone._atom2 = atom_map[self._atom2.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._atom1}, {self._carbon2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
