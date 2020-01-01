from .. import FunctionalGroup_


class Thiol(FunctionalGroup_):
    """
    Represents a thiol functional group.


    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][sulfur][hydrogen]``.

    """

    def __init__(self, sulfur, hydrogen, atom, bonders, deleters):
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_sulfur(self):
        return self._sulfur.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (self._sulfur, self._hydrogen, self._atom)
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._sulfur = atom_map[self._sulfur.id]
        clone._hydrogen = atom_map[self._hydrogen.id]
        clone._atom = atom_map[self._atom.id]
        return clone
