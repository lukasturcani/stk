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
        return self._fluorine.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (self._fluorine, self._atom)
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._fluorine = atom_map[self._fluorine.id]
        clone._atom = atom_map[self._atom.id]
        return clone
