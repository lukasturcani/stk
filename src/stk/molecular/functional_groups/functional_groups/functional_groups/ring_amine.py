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
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._hydrogen3 = hydrogen3
        self._carbon1 = carbon1
        self._carbon2 = carbon2
        self._carbon3 = carbon3
        atoms = (
            nitrogen,
            hydrogen1,
            hydrogen2,
            carbon1,
            carbon2,
            hydrogen3,
            carbon3,
        )
        super().__init__(atoms, (), ())

    def get_nitrogen(self):
        return self._nitrogen

    def get_hydrogen1(self):
        return self._hydrogen1

    def get_hydrogen2(self):
        return self._hydrogen2

    def get_carbon1(self):
        return self._carbon1

    def get_carbon2(self):
        return self._carbon2

    def get_hydrogen3(self):
        return self._hydrogen3

    def get_carbon3(self):
        return self._carbon3

    def clone(self):
        clone = super().clone()
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
        clone._hydrogen2 = self._hydrogen2
        clone._hydrogen3 = self._hydrogen3
        clone._carbon1 = self._carbon1
        clone._carbon2 = self._carbon2
        clone._carbon3 = self._carbon3
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        clone._hydrogen3 = atom_map.get(
            self._hydrogen3.get_id(),
            self._hydrogen3,
        )
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        clone._carbon3 = atom_map.get(
            self._carbon3.get_id(),
            self._carbon3,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._carbon1}, {self._carbon2}, {self._hydrogen3}, '
            f'{self._carbon3})'
        )
