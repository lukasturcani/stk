from .. import FunctionalGroup_


class Fluoro(FunctionalGroup_):
    """
    Represents a fluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine][atom]``.

    """

    def __init__(self, fluorine, atom, bonders, deleters):
        self._fluorine = fluorine
        self._atom = atom
        super().__init__((fluorine, atom), bonders, deleters)

    def get_fluorine(self):
        return self._fluorine

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._fluorine = self._fluorine
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._fluorine = atom_map.get(
            self._fluorine.get_id(),
            self._fluorine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
