from .. import FunctionalGroup_


class Bromo(FunctionalGroup_):
    """
    Represents a bromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine][atom]``.

    """

    def __init__(self, bromine, atom, bonders, deleters):
        self._bromine = bromine
        self._atom = atom
        super().__init__((bromine, atom), bonders, deleters)

    def get_bromine(self):
        return self._bromine

    def get_atom(self):
        return self._atom

    def clone(self):
        clone = super().clone()
        clone._bromine = self._bromine
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._bromine = atom_map.get(
            self._bromine.get_id(),
            self._bromine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
