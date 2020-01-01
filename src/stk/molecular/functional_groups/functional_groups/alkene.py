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
        self._carbon1 = carbon1
        self._atom1 = atom1
        self._atom2 = atom2
        self._carbon2 = carbon2
        self._atom3 = atom3
        self._atom4 = atom4
        atoms = (carbon1, atom1, atom2, carbon2, atom3, atom4)
        super().__init__(atoms, bonders, deleters)

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
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon1 = atom_map[self._carbon1.id]
        clone._atom1 = atom_map[self._atom1.id]
        clone._atom2 = atom_map[self._atom2.id]
        clone._carbon2 = atom_map[self._carbon2.id]
        clone._atom3 = atom_map[self._atom3.id]
        clone._atom4 = atom_map[self._atom4.id]
        return clone
