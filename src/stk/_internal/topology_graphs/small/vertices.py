import logging

import numpy as np

from stk._internal.building_block import BuildingBlock
from stk._internal.topology_graphs.vertex import Vertex
from stk._internal.utilities.utilities import get_acute_vector

from ..edge import Edge
from ..utilities import _EdgeSorter, _FunctionalGroupSorter

logger = logging.getLogger(__name__)


class CoreVertex(Vertex):
    """
    Represents a vertex in the core of an ncore topology graph.

    """

    def place_building_block(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> np.ndarray:
        # Sets building block in XY plane, then aligns 0th functional
        # group with the 0th edge.
        assert building_block.get_num_functional_groups() > 2, (
            f"{building_block} needs to have more than 2 functional "
            "groups but has "
            f"{building_block.get_num_functional_groups()}."
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = get_acute_vector(
            reference=core_centroid - self._position,
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=np.array([0.0, 0.0, 1.0]),
            origin=self._position,
        )
        (fg,) = building_block.get_functional_groups(0)
        fg_centroid = building_block.get_centroid(fg.get_placer_ids())
        edge_position = edges[0].get_position()
        return building_block.with_rotation_to_minimize_angle(
            start=fg_centroid - self._position,
            target=edge_position - self._position,
            axis=np.array([0, 0, 1], dtype=np.float64),
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> dict[int, int]:
        # The idea is to order the functional groups in building_block
        # by their angle with the vector running from the placer
        # centroid to fg0, going in the clockwise direction.
        # The edges are also ordered by their angle with the vector
        # running from the edge centroid to the aligner_edge,
        # going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg_sorter = _FunctionalGroupSorter(building_block)
        edge_sorter = _EdgeSorter(
            edges=edges,
            aligner_edge=edges[0],
            axis=fg_sorter.get_axis(),
        )
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(
                fg_sorter.get_items(),
                edge_sorter.get_items(),
            )
        }


class SubstituentVertex(Vertex):
    """
    Represents a vertex to be bound to core.

    """

    def place_building_block(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> np.ndarray:
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
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(
        self,
        building_block: BuildingBlock,
        edges: tuple[Edge, ...],
    ) -> dict[int, int]:
        return {0: edges[0].get_id()}
