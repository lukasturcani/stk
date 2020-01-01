from .. import FunctionalGroup_


class SecondaryAmine(FunctionalGroup_):
    """
    Represents a secondary amine functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom1][nitrogen]([hydrogen])[atom2]``.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen,
        atom1,
        atom2,
        bonders,
        deleters,
    ):
        self._nitrogen = nitrogen
        self._hydrogen = hydrogen
        self._atom1 = atom1
        self._atom2 = atom2
        atoms = (nitrogen, hydrogen, atom1, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom1(self):
        return self._atom1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._nitrogen, self._hydrogen, self._atom1, self._atom2
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._nitrogen = atom_map[self._nitrogen.id]
        clone._hydrogen = atom_map[self._hydrogen.id]
        clone._atom1 = atom_map[self._atom1.id]
        clone._atom2 = atom_map[self._atom2.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen}, {self._atom1}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
