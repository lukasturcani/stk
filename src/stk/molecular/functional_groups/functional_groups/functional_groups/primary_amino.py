from .. import FunctionalGroup_


class PrimaryAmino(FunctionalGroup_):
    """
    Represents a primary amino functional group.

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
        return self._nitrogen

    def get_hydrogen1(self):
        return self._hydrogen1

    def get_hydrogen2(self):
        return self._hydrogen2

    def get_atom(self):
        return self._atom

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters}'
            ')'
        )

    def clone(self):
        clone = super().clone()
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone
