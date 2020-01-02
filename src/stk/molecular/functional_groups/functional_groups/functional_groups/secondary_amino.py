from .. import FunctionalGroup_


class SecondaryAmino(FunctionalGroup_):
    """
    Represents a secondary amino functional group.

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
        return self._nitrogen

    def get_hydrogen(self):
        return self._hydrogen

    def get_atom1(self):
        return self._atom1

    def get_atom2(self):
        return self._atom2

    def clone(self):
        clone = super().clone()
        clone._nitrogen = self._nitrogen
        clone._hydrogen = self._hydrogen
        clone._atom1 = self._atom1
        clone._atom2 = self._atom2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen}, {self._atom1}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
