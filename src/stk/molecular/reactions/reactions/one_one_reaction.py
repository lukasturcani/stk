"""
One-One Reaction
================

"""

from ...bonds import Bond
from .reaction import Reaction


class OneOneReaction(Reaction):
    """
    A reaction between two functional groups, each with 1 bonder atom.

    The reaction creates a bond between the *bonder* atoms, and deletes
    any *deleter* atoms.

    """

    def __init__(
        self,
        functional_group1,
        functional_group2,
        bond_order,
        periodicity,
    ):
        """
        Initialize a :class:`.OneOneReaction` instance.

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
        bonder1 = next(self._functional_group1.get_bonders())
        bonder2 = next(self._functional_group2.get_bonders())
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
