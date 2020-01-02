from .. import FunctionalGroup_


class Amide(FunctionalGroup_):
    """
    Represents an amide functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        nitrogen,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon

    def get_oxygen(self):
        return self._oxygen

    def get_nitrogen(self):
        return self._nitrogen

    def get_hydrogen1(self):
        return self._hydrogen1

    def get_hydrogen2(self):
        return self._hydrogen2

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._oxygen = atom_map.get(
            self._oxygen.get_id(),
            self._oxygen,
        )
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

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._nitrogen}, '
            f'{self._hydrogen1}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
