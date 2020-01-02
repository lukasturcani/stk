from ..reaction import Reaction
from ...bond import Bond


class OneOneReaction(Reaction):
    """
    A reaction between two functional groups, each with 1 bonder atom.

    """

    def __init__(
        self,
        functional_group1,
        functional_group2,
        bond_order,
        periodicity,
    ):
        self._functional_group1 = functional_group1.clone()
        self._functional_group2 = functional_group2.clone()
        self._bond_order = bond_order
        self._periodicity = periodicity

    def get_new_atoms(self):
        return
        yield

    def get_new_bonds(self):
        bonder1 = next(self._functional_group1.get_bonders())
        bonder2 = next(self._functional_group2.get_bonders())
        yield Bond(
            atom1=bonder1,
            atom2=bonder2,
            order=self._bond_order,
            periodicity=self._periodicity,
        )

    def get_deleted_atoms(self):
        yield from self._functional_group1.get_deleters()
        yield from self._functional_group2.get_deleters()
