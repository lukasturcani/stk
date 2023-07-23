"""
Polymer and Small Molecule Vertices
===================================

"""

import logging

import numpy as np

from stk._internal.topology_graphs.vertex import Vertex
from stk._internal.utilities.utilities import get_acute_vector

from ..utilities import _EdgeSorter, _FunctionalGroupSorter

logger = logging.getLogger(__name__)


class LinearVertex(Vertex):
    """
    Represents a vertex in the middle of a linear polymer chain.

    """

    def __init__(self, id, position, flip):
        """
        Initialize a :class:`.LinearVertex` instance.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`numpy.ndarray`
            The position of the vertex.

        flip : :class:`bool`
            If ``True`` any building block placed by the vertex will
            have its orientation along the chain flipped.

        """

        super().__init__(id, position)
        self._flip = flip

    def get_flip(self):
        """
        Return ``True`` if the vertex flips building blocks it places.

        Returns
        -------
        :class:`bool`
            ``True`` if the vertex flips building blocks it places.

        """

        return self._flip

    def clone(self):
        clone = super().clone()
        clone._flip = self._flip
        return clone

    def place_building_block(self, building_block, edges):
        assert building_block.get_num_functional_groups() == 2, (
            f"{building_block} needs to have exactly 2 functional "
            "groups but has "
            f"{building_block.get_num_functional_groups()}."
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
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
        fg1_id, fg2_id = self._sort_functional_groups(building_block)
        edge1_id, edge2_id = self._sort_edges(edges)
        return {
            fg1_id: edge1_id,
            fg2_id: edge2_id,
        }

    @staticmethod
    def _sort_functional_groups(building_block):
        fg1, fg2 = building_block.get_functional_groups()
        x1, y1, z1 = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        x2, y2, z2 = building_block.get_centroid(
            atom_ids=fg2.get_placer_ids(),
        )
        return (0, 1) if x1 < x2 else (1, 0)

    @staticmethod
    def _sort_edges(edges):
        edge1, edge2 = edges
        x1, y1, z1 = edge1.get_position()
        x2, y2, z2 = edge2.get_position()
        if x1 < x2:
            return edge1.get_id(), edge2.get_id()
        else:
            return edge2.get_id(), edge1.get_id()

    def __str__(self):
        return (
            f"Vertex(id={self._id}, "
            f"position={self._position.tolist()}, "
            f"flip={self._flip})"
        )


class TerminalVertex(LinearVertex):
    """
    Represents a vertex at the end of a polymer chain.

    Do not instantiate this class directly, use :class:`.HeadVertex`
    or :class:`.TailVertex` instead.

    """

    def place_building_block(self, building_block, edges):
        if (
            building_block.get_num_functional_groups() != 1
            and building_block.get_num_placers() > 1
        ):
            return super().place_building_block(building_block, edges)

        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg, *_ = building_block.get_functional_groups()
        fg_centroid = building_block.get_centroid(
            atom_ids=fg.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        return building_block.with_rotation_between_vectors(
            start=fg_centroid - core_centroid,
            # _cap_direction is defined by a subclass.
            target=[self._cap_direction, 0, 0],
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        if building_block.get_num_functional_groups() == 2:
            functional_groups = self._sort_functional_groups(
                building_block=building_block,
            )
            index = 1 if self._cap_direction == 1 else 0
            return {functional_groups[index]: edges[0].get_id()}

        elif building_block.get_num_functional_groups() == 1:
            return {0: edges[0].get_id()}

        else:
            raise ValueError(
                "The building block of a polymer "
                "must have 1 or 2 functional groups."
            )


class HeadVertex(TerminalVertex):
    """
    Represents a vertex at the head of a polymer chain.

    """

    # The direction to use if the building block placed on the
    # vertex only has 1 FunctionalGroup.
    _cap_direction = 1


class TailVertex(TerminalVertex):
    """
    Represents a vertex at the tail of a polymer chain.

    """

    # The direction to use if the building block placed on the
    # vertex only has 1 FunctionalGroup.
    _cap_direction = -1


class UnaligningVertex(LinearVertex):
    """
    Just places a building block, does not align.

    """

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {fg_id: edge.get_id() for fg_id, edge in enumerate(edges)}


class CoreVertex(Vertex):
    """
    Represents a vertex in the core of an ncore topology graph.

    """

    def place_building_block(self, building_block, edges):
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
            target=[0, 0, 1],
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

    # def map_functional_groups_to_edges(self, building_block, edges):
    #     return {fg_id: edge.get_id() for fg_id, edge in enumerate(edges)}

    def map_functional_groups_to_edges(self, building_block, edges):
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
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        return {0: edges[0].get_id()}
