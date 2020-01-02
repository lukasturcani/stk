from .. import FunctionalGroup_


class Alcohol(FunctionalGroup_):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(self, oxygen, hydrogen, atom, bonders, deleters):
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (oxygen, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_oxygen(self):
        return self._oxygen

    def get_hydrogen(self):
        return self._hydrogen

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
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

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._oxygen}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
