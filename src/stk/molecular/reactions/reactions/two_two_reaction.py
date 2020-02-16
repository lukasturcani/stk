import numpy as np
import itertools as it
from scipy.spatial.distance import euclidean

from .reaction import Reaction
from ...bond import Bond


class TwoTwoReaction(Reaction):
    """
    A reaction between two functional groups, each with 2 bonder atoms.

    """

    def __init__(
        self,
        position_matrix,
        functional_group1,
        functional_group2,
        bond_order,
        periodicity,
    ):
        self._position_matrix = np.array(position_matrix)
        self._functional_group1 = functional_group1.clone()
        self._functional_group2 = functional_group2.clone()
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
