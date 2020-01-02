from .. import FunctionalGroup_


class Iodo(FunctionalGroup_):
    """
    Represents an iodo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[iodine][atom]``.

    """

    def __init__(self, iodine, atom, bonders, deleters):
        self._iodine = iodine
        self._atom = atom
        super().__init__((iodine, atom), bonders, deleters)

    def get_iodine(self):
        return self._iodine

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._iodine = self._iodine
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._iodine = atom_map.get(
            self._iodine.get_id(),
            self._iodine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._iodine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
