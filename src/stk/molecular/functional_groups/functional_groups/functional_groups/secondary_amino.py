from .. import FunctionalGroup_


class SecondaryAmino(FunctionalGroup_):
    """
    Represents a secondary amino functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom1][nitrogen]([hydrogen])[atom2]``.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen,
        atom1,
        atom2,
        bonders,
        deleters,
    ):
        atom_map = {
            nitrogen.get_id(): nitrogen.clone(),
            hydrogen.get_id(): hydrogen.clone(),
            atom1.get_id(): atom1.clone(),
            atom2.get_id(): atom2.clone(),
        }
        self._nitrogen = atom_map[nitrogen.get_id()]
        self._hydrogen = atom_map[hydrogen.get_id()]
        self._atom1 = atom_map[atom1.get_id()]
        self._atom2 = atom_map[atom2.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
        )

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom1(self):
        return self._atom1.clone()

    def get_atom2(self):
        return self._atom2.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._nitrogen, self._hydrogen, self._atom1, self._atom2
        )
        for atom in atoms:
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._nitrogen = atom_map[self._nitrogen.get_id()]
        clone._hydrogen = atom_map[self._hydrogen.get_id()]
        clone._atom1 = atom_map[self._atom1.get_id()]
        clone._atom2 = atom_map[self._atom2.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen}, {self._atom1}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
