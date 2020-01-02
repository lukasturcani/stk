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
            fluorine1.get_id(): fluorine1.clone(),
            atom1.get_id(): atom1.clone(),
            fluorine2.get_id(): fluorine2.clone(),
            atom2.get_id(): atom2.clone(),
        }
        self._atom1 = atom_map[atom1.get_id()]
        self._fluorine1 = atom_map[fluorine1.get_id()]
        self._atom2 = atom_map[atom2.get_id()]
        self._fluorine2 = atom_map[fluorine2.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
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
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._atom1 = atom_map[self._atom1.get_id()]
        clone._fluorine1 = atom_map[self._fluorine1.get_id()]
        clone._atom2 = atom_map[self._atom2.get_id()]
        clone._fluorine2 = atom_map[self._fluorine2.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine1}, {self._atom1}, {self._fluorine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
