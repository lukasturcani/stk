"""
Functional Group Sorter
=======================

"""

import numpy as np

from stk.utilities import get_acute_vector
from ...building_block import BuildingBlock
from .angle_sorter import AngleSorter


__all__ = (
    'FunctionalGroupSorter',
)


class FunctionalGroupSorter(AngleSorter[int]):
    """
    Sorts functional groups according to their angle.

    """

    __slots__ = [
        '_items',
        '_reference',
        '_axis',
        '_placer_centroid',
        '_building_block',
    ]

    def __init__(
        self,
        building_block: BuildingBlock,
    ) -> None:
        """
        Initialize a :class:`.FunctionalGroupSorter` instance.

        Parameters:

            building_block:
                The building block, whose functional groups are to be
                sorted.

        """

        self._building_block = building_block
        fg0_position = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        self._placer_centroid = placer_centroid = (
            building_block.get_centroid(
                atom_ids=building_block.get_placer_ids(),
            )
        )
        fg0_direction = fg0_position - placer_centroid
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        axis = np.cross(
            fg0_direction,
            get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(),
            ),
        )
        axis.setflags(write=False)
        super().__init__(
            items=range(building_block.get_num_functional_groups()),
            reference=fg0_direction,
            axis=axis,
        )

    def _get_vector(
        self,
        item: int,
    ) -> np.ndarray:
        building_block = self._building_block
        fg, = building_block.get_functional_groups(item)
        fg_position = building_block.get_centroid(fg.get_placer_ids())
        return fg_position - self._placer_centroid
