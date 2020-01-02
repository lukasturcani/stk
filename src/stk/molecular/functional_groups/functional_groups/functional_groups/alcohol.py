from .. import FunctionalGroup_


class Alcohol(FunctionalGroup_):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(self, oxygen, hydrogen, atom, bonders, deleters):
        atom_map = {
            oxygen.get_id(): oxygen.clone(),
            hydrogen.get_id(): hydrogen.clone(),
            atom.get_id(): atom.clone(),
        }
        self._oxygen = atom_map[oxygen.get_id()]
        self._hydrogen = atom_map[hydrogen.get_id()]
        self._atom = atom_map[atom.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
        )

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (self._oxygen, self._hydrogen, self._atom)
        for atom in atoms:
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._oxygen = atom_map[self._oxygen.get_id()]
        clone._hydrogen = atom_map[self._hydrogen.get_id()]
        clone._atom = atom_map[self._atom.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._oxygen}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
