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
        self._atom1 = atom1
        self._fluorine1 = fluorine1
        self._atom2 = atom2
        self._fluorine2 = fluorine2
        atoms = (atom1, fluorine1, atom2, fluorine2)
        super().__init__(atoms, bonders, deleters)

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
