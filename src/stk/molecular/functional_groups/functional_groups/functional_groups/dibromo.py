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
            bromine1.get_id(): bromine1.clone(),
            atom1.get_id(): atom1.clone(),
            bromine2.get_id(): bromine2.clone(),
            atom2.get_id(): atom2.clone(),
        }
        self._atom1 = atom_map[atom1.get_id()]
        self._bromine1 = atom_map[bromine1.get_id()]
        self._atom2 = atom_map[atom2.get_id()]
        self._bromine2 = atom_map[bromine2.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
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
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._atom1 = atom_map[self._atom1.get_id()]
        clone._bromine1 = atom_map[self._bromine1.get_id()]
        clone._atom2 = atom_map[self._atom2.get_id()]
        clone._bromine2 = atom_map[self._bromine2.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine1}, {self._atom1}, {self._bromine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
