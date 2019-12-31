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
        super().__init__(
            atoms=(
                nitrogen,
                hydrogen1,
                hydrogen2,
                carbon1,
                carbon2,
                hydrogen3,
                carbon3,
            ),
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
