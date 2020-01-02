from .. import FunctionalGroup_


class CarboxylicAcid(FunctionalGroup_):
    """
    Represents a carboxylic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen1])[oxygen2][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen1,
        oxygen2,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        atom_map = {
            carbon.id: carbon.clone(),
            oxygen1.id: oxygen1.clone(),
            oxygen2.id: oxygen2.clone(),
            hydrogen.id: hydrogen.clone(),
            atom.id: atom.clone(),
        }
        self._carbon = atom_map[carbon.id]
        self._oxygen1 = atom_map[oxygen1.id]
        self._oxygen2 = atom_map[oxygen2.id]
        self._hydrogen = atom_map[hydrogen.id]
        self._atom = atom_map[atom.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.id] for a in bonders),
            deleters=tuple(atom_map[a.id] for a in deleters),
        )

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen1(self):
        return self._oxygen1.clone()

    def get_oxygen2(self):
        return self._oxygen2.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._carbon,
            self._oxygen1,
            self._oxygen2,
            self._hydrogen,
            self._atom,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.id]
        clone._oxygen1 = atom_map[self._oxygen1.id]
        clone._oxygen2 = atom_map[self._oxygen2.id]
        clone._hydrogen = atom_map[self._hydrogen.id]
        clone._atom = atom_map[self._atom.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen1}, {self._oxygen2}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
