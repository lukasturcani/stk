from .. import FunctionalGroup_


class Thioacid(FunctionalGroup_):
    """
    Represents a thioacid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[sulfur][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        sulfur,
        hydrogen,
        atom,
        bonders,
        deleters
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon

    def get_oxygen(self):
        return self._oxygen

    def get_sulfur(self):
        return self._sulfur

    def get_hydrogen(self):
        return self._hydrogen

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._sulfur = self._sulfur
        clone._hydrogen = self._hydrogen
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
        clone._sulfur = atom_map.get(
            self._sulfur.get_id(),
            self._sulfur,
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

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._sulfur}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
