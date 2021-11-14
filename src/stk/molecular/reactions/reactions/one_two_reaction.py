"""
One-Two Reaction
================

"""

import itertools as it

from ...bonds import Bond
from .reaction import Reaction


class OneTwoReaction(Reaction):
    """
    A reaction between two functional groups.

    One functional group has 1 *bonder* atom, A, and the other has 2, B
    and C. Two bonds are created, one between A and B and one between
    A and C. Any deleter atoms are removed.

    """

    def __init__(
        self,
        functional_group1,
        functional_group2,
        bond_order,
        periodicity,
    ):
        """
        Initialize a :class:`.OneTwoReaction` instance.

        Parameters
        ----------
        functional_group1 : :class:`.GenericFunctionalGroup`
            The first functional group in the reaction.

        functional_group2 : :class:`.GenericFunctionalGroup`
            The second functional group in the reaction.

        bond_order : :class:`int`
            The bond order of the bond created by the reaction.

        periodicity : :class:`tuple` of :class:`int`
            The periodicity of the bond created by the reaction.

        """

        self._functional_group1 = functional_group1
        self._functional_group2 = functional_group2
        self._bond_order = bond_order
        self._periodicity = periodicity

    def _get_new_atoms(self):
        return
        yield

    def _get_new_bonds(self):
        bonders1 = self._functional_group1.get_bonders()
        bonders2 = self._functional_group2.get_bonders()
        for bonder1, bonder2 in it.product(bonders1, bonders2):
            yield Bond(
                atom1=bonder1,
                atom2=bonder2,
                order=self._bond_order,
                periodicity=self._periodicity,
            )

    def _get_deleted_atoms(self):
        yield from self._functional_group1.get_deleters()
        yield from self._functional_group2.get_deleters()

    def _get_deleted_bonds(self):
        return
        yield
