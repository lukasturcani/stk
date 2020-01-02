from .. import FunctionalGroup_


class Iodo(FunctionalGroup_):
    """
    Represents an iodo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[iodine][atom]``.

    """

    def __init__(self, iodine, atom, bonders, deleters):
        atom_map = {
            iodine.get_id(): iodine.clone(),
            atom.get_id(): atom.clone(),
        }
        self._iodine = atom_map[iodine.get_id()]
        self._atom = atom_map[atom.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
        )

    def get_iodine(self):
        return self._iodine.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (self._iodine, self._atom)
        for atom in atoms:
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._iodine = atom_map[self._iodine.get_id()]
        clone._atom = atom_map[self._atom.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._iodine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
