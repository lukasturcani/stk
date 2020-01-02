from .. import FunctionalGroup_


class BoronicAcid(FunctionalGroup_):
    """
    Represents a boronic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][boron]([oxygen1][hydrogen1])[oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        boron,
        oxygen1,
        hydrogen1,
        oxygen2,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        atom_map = {
            boron.id: boron.clone(),
            oxygen1.id: oxygen1.clone(),
            hydrogen1.id: hydrogen1.clone(),
            oxygen2.id: oxygen2.clone(),
            hydrogen2.id: hydrogen2.clone(),
            atom.id: atom.clone(),
        }
        self._boron = atom_map[boron.id]
        self._oxygen1 = atom_map[oxygen1.id]
        self._hydrogen1 = atom_map[hydrogen1.id]
        self._oxygen2 = atom_map[oxygen2.id]
        self._hydrogen2 = atom_map[hydrogen2.id]
        self._atom = atom_map[atom.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.id] for a in bonders),
            deleters=tuple(atom_map[a.id] for a in deleters),
        )

    def get_boron(self):
        return self._boron.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

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
            self._boron,
            self._oxygen1,
            self._hydrogen1,
            self._oxygen2,
            self._hydrogen2,
            self._atom,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._boron = atom_map[self._boron.id]
        clone._oxygen1 = atom_map[self._oxygen1.id]
        clone._hydrogen1 = atom_map[self._hydrogen1.id]
        clone._oxygen2 = atom_map[self._oxygen2.id]
        clone._hydrogen2 = atom_map[self._hydrogen2.id]
        clone._atom = atom_map[self._atom.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._boron}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._oxygen2}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
