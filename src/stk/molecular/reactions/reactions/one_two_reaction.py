import itertools as it

from ..reaction import Reaction
from ...bond import Bond


class OneTwoReaction(Reaction):
    """
    A reaction between two functional groups.

    One functional group has 1 bonder atom and the other has 2.

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
        bonders1 = self._functional_group1.get_bonders()
        bonders2 = self._functional_group2.get_bonders()
        for bonder1, bonder2 in it.product(bonders1, bonders2):
            yield Bond(
                atom1=bonder1,
                atom2=bonder2,
                order=self._bond_order,
                periodicity=self._periodicity,
            )

    def get_deleted_atoms(self):
        yield from self._functional_group1.get_deleters()
        yield from self._functional_group2.get_deleters()
