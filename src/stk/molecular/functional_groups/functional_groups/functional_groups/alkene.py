from .. import FunctionalGroup_


class Alkene(FunctionalGroup_):
    """
    Represents an alkene functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon1]([atom1])([atom2])=[carbon2]([atom3])[atom4]``.

    """

    def __init__(
        self,
        carbon1,
        atom1,
        atom2,
        carbon2,
        atom3,
        atom4,
        bonders,
        deleters,
    ):
        atom_map = {
            carbon1.get_id(): carbon1.clone(),
            atom1.get_id(): atom1.clone(),
            atom2.get_id(): atom2.clone(),
            carbon2.get_id(): carbon2.clone(),
            atom3.get_id(): atom3.clone(),
            atom4.get_id(): atom4.clone(),
        }
        self._carbon1 = atom_map[carbon1.get_id()]
        self._atom1 = atom_map[atom1.get_id()]
        self._atom2 = atom_map[atom2.get_id()]
        self._carbon2 = atom_map[carbon2.get_id()]
        self._atom3 = atom_map[atom3.get_id()]
        self._atom4 = atom_map[atom4.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
        )

    def get_carbon1(self):
        return self._carbon1.clone()

    def get_atom1(self):
        return self.atom1.clone()

    def get_atom2(self):
        return self.atom2.clone()

    def get_carbon2(self):
        return self._carbon2.clone()

    def get_atom3(self):
        return self._atom3.clone()

    def get_atom4(self):
        return self._atom4.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._carbon1,
            self._atom1,
            self._atom2,
            self._carbon2,
            self._atom3,
            self._atom4,
        )
        for atom in atoms:
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon1 = atom_map[self._carbon1.get_id()]
        clone._atom1 = atom_map[self._atom1.get_id()]
        clone._atom2 = atom_map[self._atom2.get_id()]
        clone._carbon2 = atom_map[self._carbon2.get_id()]
        clone._atom3 = atom_map[self._atom3.get_id()]
        clone._atom4 = atom_map[self._atom4.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._atom1}, {self._atom2}, '
            f'{self._carbon2}, {self._atom3}, {self._atom4}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
