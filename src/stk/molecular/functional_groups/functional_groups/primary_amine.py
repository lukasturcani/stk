from .. import FunctionalGroup_


class PrimaryAmine(FunctionalGroup_):
    """
    Represents a primary amine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (nitrogen, hydrogen1, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters}'
            ')'
        )

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._nitrogen,
            self._hydrogen1,
            self._hydrogen2,
            self._atom,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._nitrogen = atom_map[self._nitrogen.id]
        clone._hydrogen1 = atom_map[self._hydrogen1.id]
        clone._hydrogen2 = atom_map[self._hydrogen2.id]
        clone._atom = atom_map[self._atom.id]
        return clone
