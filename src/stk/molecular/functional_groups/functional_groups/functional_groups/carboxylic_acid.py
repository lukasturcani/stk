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
            carbon.get_id(): carbon.clone(),
            oxygen1.get_id(): oxygen1.clone(),
            oxygen2.get_id(): oxygen2.clone(),
            hydrogen.get_id(): hydrogen.clone(),
            atom.get_id(): atom.clone(),
        }
        self._carbon = atom_map[carbon.get_id()]
        self._oxygen1 = atom_map[oxygen1.get_id()]
        self._oxygen2 = atom_map[oxygen2.get_id()]
        self._hydrogen = atom_map[hydrogen.get_id()]
        self._atom = atom_map[atom.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
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
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.get_id()]
        clone._oxygen1 = atom_map[self._oxygen1.get_id()]
        clone._oxygen2 = atom_map[self._oxygen2.get_id()]
        clone._hydrogen = atom_map[self._hydrogen.get_id()]
        clone._atom = atom_map[self._atom.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen1}, {self._oxygen2}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
