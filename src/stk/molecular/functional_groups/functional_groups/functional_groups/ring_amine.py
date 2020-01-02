from .. import FunctionalGroup_


class RingAmine(FunctionalGroup_):
    """
    Represents an amine bonded to a ring.

    The structure of the functional group is given by the pseudo-SMILES
    ``[hydrogen1][nitrogen]([hydrogen2])[carbon1][carbon2]
    ([hydrogen3])[carbon3]``.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen1,
        hydrogen2,
        carbon1,
        carbon2,
        hydrogen3,
        carbon3,
    ):
        atom_map = {
            nitrogen.id: nitrogen.clone(),
            hydrogen1.id: hydrogen1.clone(),
            hydrogen2.id: hydrogen2.clone(),
            carbon1.id: carbon1.clone(),
            carbon2.id: carbon2.clone(),
            hydrogen3.id: hydrogen3.clone(),
            carbon3.id: carbon3.clone(),
        }
        self._nitrogen = atom_map[nitrogen.id]
        self._hydrogen1 = atom_map[hydrogen1.id]
        self._hydrogen2 = atom_map[hydrogen2.id]
        self._hydrogen3 = atom_map[hydrogen3.id]
        self._carbon1 = atom_map[carbon1.id]
        self._carbon2 = atom_map[carbon2.id]
        self._carbon3 = atom_map[carbon3.id]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=(),
            deleters=(),
        )

    def get_nitrogen(self):
        return self._nitrogen.clone()

    def get_hydrogen1(self):
        return self._hydrogen1.clone()

    def get_hydrogen2(self):
        return self._hydrogen2.clone()

    def get_carbon1(self):
        return self._carbon1.clone()

    def get_carbon2(self):
        return self._carbon2.clone()

    def get_hydrogen3(self):
        return self._hydrogen3.clone()

    def get_carbon3(self):
        return self._carbon3.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._nitrogen,
            self._hydrogen1,
            self._hydrogen2,
            self._hydrogen3,
            self._carbon1,
            self._carbon2,
            self._carbon3,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._nitrogen = atom_map[self._nitrogen.id]
        clone._hydrogen1 = atom_map[self._hydrogen1.id]
        clone._hydrogen2 = atom_map[self._hydrogen2.id]
        clone._hydrogen3 = atom_map[self._hydrogen3.id]
        clone._carbon1 = atom_map[self._carbon1.id]
        clone._carbon2 = atom_map[self._carbon2.id]
        clone._carbon3 = atom_map[self._carbon3.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._carbon1}, {self._carbon2}, {self._hydrogen3}, '
            f'{self._carbon3})'
        )
