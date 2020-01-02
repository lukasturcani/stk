from .. import FunctionalGroup_


class Amide(FunctionalGroup_):
    """
    Represents an amide functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        nitrogen,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        atom_map = {
            carbon.get_id(): carbon.clone(),
            oxygen.get_id(): oxygen.clone(),
            nitrogen.get_id(): nitrogen.clone(),
            hydrogen1.get_id(): hydrogen1.clone(),
            hydrogen2.get_id(): hydrogen2.clone(),
            atom.get_id(): atom.clone(),
        }
        self._carbon = atom_map[carbon.get_id()]
        self._oxygen = atom_map[oxygen.get_id()]
        self._nitrogen = atom_map[nitrogen.get_id()]
        self._hydrogen1 = atom_map[hydrogen1.get_id()]
        self._hydrogen2 = atom_map[hydrogen2.get_id()]
        self._atom = atom_map[atom.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
        )

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._carbon,
            self._oxygen,
            self._nitrogen,
            self._hydrogen1,
            self._hydrogen2,
            self._atom,
        )
        for atom in atoms:
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.get_id()]
        clone._oxygen = atom_map[self._oxygen.get_id()]
        clone._nitrogen = atom_map[self._nitrogen.get_id()]
        clone._hydrogen1 = atom_map[self._hydrogen1.get_id()]
        clone._hydrogen2 = atom_map[self._hydrogen2.get_id()]
        clone._atom = atom_map[self._atom.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._nitrogen}, '
            f'{self._hydrogen1}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
