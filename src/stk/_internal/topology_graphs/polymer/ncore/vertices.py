"""
Linear Polymer Vertices
=======================

"""

import logging

from stk._internal.topology_graphs.vertex import Vertex
from ...utilities import _EdgeSorter, _FunctionalGroupSorter

logger = logging.getLogger(__name__)


class CoreVertex(Vertex):
    """
    Represents a vertex in the middle of a linear polymer chain.

    """

    def place_building_block(self, building_block, edges):
        assert building_block.get_num_functional_groups() > 2, (
            f"{building_block} needs to have more than 2 functional "
            "groups but has "
            f"{building_block.get_num_functional_groups()}."
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        return building_block.get_position_matrix()
        fg1, fg2 = building_block.get_functional_groups()
        fg1_position = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        fg2_position = building_block.get_centroid(
            atom_ids=fg2.get_placer_ids(),
        )
        return building_block.with_rotation_between_vectors(
            start=fg2_position - fg1_position,
            target=[-1 if self._flip else 1, 0, 0],
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {fg_id: edge.get_id() for fg_id, edge in enumerate(edges)}

    # def map_functional_groups_to_edges(self, building_block, edges):
    #     # The idea is to order the functional groups in building_block
    #     # by their angle with the vector running from the placer
    #     # centroid to fg0, going in the clockwise direction.
    #     # The edges are also ordered by their angle with the vector
    #     # running from the edge centroid to the aligner_edge,
    #     # going in the clockwise direction.
    #     #
    #     # Once the fgs and edges are ordered, zip and assign them.

    #     fg_sorter = _FunctionalGroupSorter(building_block)
    #     edge_sorter = _EdgeSorter(
    #         edges=edges,
    #         aligner_edge=edges[self._aligner_edge],
    #         axis=fg_sorter.get_axis(),
    #     )
    #     return {
    #         fg_id: edge.get_id()
    #         for fg_id, edge in zip(
    #             fg_sorter.get_items(),
    #             edge_sorter.get_items(),
    #         )
    #     }


class TerminalVertex(Vertex):
    """
    Represents a vertex to be bound to core.

    """

    def place_building_block(self, building_block, edges):
        assert building_block.get_num_functional_groups() == 1, (
            f"{building_block} needs to have 1 functional "
            "groups but has "
            f"{building_block.get_num_functional_groups()}."
        )

        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_core_atom_ids(),
        )

        (fg,) = building_block.get_functional_groups()
        fg_centroid = building_block.get_centroid(
            atom_ids=fg.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        edge_centroid = sum(edge.get_position() for edge in edges) / len(edges)
        return building_block.with_rotation_between_vectors(
            start=(fg_centroid - core_centroid),
            # _cap_direction is defined by a subclass.
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {0: edges[0].get_id()}
