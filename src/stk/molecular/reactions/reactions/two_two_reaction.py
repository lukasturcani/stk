"""
Two-Two Reaction
================

"""

import itertools as it

from scipy.spatial.distance import euclidean

from ...bonds import Bond
from .reaction import Reaction


class TwoTwoReaction(Reaction):
    """
    A reaction between two functional groups, each with 2 bonder atoms.

    The reaction creates the two shortest possible bonds between the
    *bonder* atoms of the two functional groups, and deletes any
    *deleter* atoms.

    """

    def __init__(
        self,
        construction_state,
        functional_group1,
        functional_group2,
        bond_order,
        periodicity,
    ):
        """
        Initialize a :class:`.TwoTwoReaction` instance.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The construction state of the molecule being constructed.

        functional_group1 : :class:`.GenericFunctionalGroup`
            The first functional group in the reaction.

        functional_group2 : :class:`.GenericFunctionalGroup`
            The second functional group in the reaction.

        bond_order : :class:`int`
            The bond order of the bonds created by the reaction.

        periodicity : :class:`tuple` of :class:`int`
            The periodicity of the bonds created by the reaction.

        """

        self._position_matrix = (
            construction_state.get_position_matrix()
        )
        self._functional_group1 = functional_group1
        self._functional_group2 = functional_group2
        self._bond_order = bond_order
        self._periodicity = periodicity

    def _get_new_atoms(self):
        return
        yield

    def _get_new_bonds(self):
        for bonder1, bonder2 in self._get_bonder_pairs():
            yield Bond(
                atom1=bonder1,
                atom2=bonder2,
                order=self._bond_order,
                periodicity=self._periodicity,
            )

    def _get_bonder_pairs(self):
        pairs = it.product(
            self._functional_group1.get_bonders(),
            self._functional_group2.get_bonders(),
        )
        sorted_pairs = sorted(pairs, key=self._pair_distance)
        bonded = set()
        for bonder1, bonder2 in sorted_pairs:
            if (
                bonder1.get_id() not in bonded
                and bonder2.get_id() not in bonded
            ):
                bonded.add(bonder1.get_id())
                bonded.add(bonder2.get_id())
                yield bonder1, bonder2

    def _pair_distance(self, bonders):
        bonder1, bonder2 = bonders
        return euclidean(
            self._position_matrix[bonder1.get_id()],
            self._position_matrix[bonder2.get_id()],
        )

    def _get_deleted_atoms(self):
        yield from self._functional_group1.get_deleters()
        yield from self._functional_group2.get_deleters()

    def _get_deleted_bonds(self):
        return
        yield
