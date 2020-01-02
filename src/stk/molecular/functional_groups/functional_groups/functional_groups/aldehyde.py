from .. import FunctionalGroup_


class Aldehyde(FunctionalGroup_):
    """
    Represents an aldehyde functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon

    def get_oxygen(self):
        return self._oxygen

    def get_hydrogen(self):
        return self._hydrogen

    def get_atom(self):
        return self._atom

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
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._hydrogen}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters})'
        )
