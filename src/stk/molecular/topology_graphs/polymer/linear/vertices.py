"""
Linear Polymer Vertices
=======================

"""

import logging

from ...topology_graph import Vertex

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
        assert (
            building_block.get_num_functional_groups() == 2
        ), (
            f'{building_block} needs to have exactly 2 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
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
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'flip={self._flip})'
        )


class TerminalVertex(LinearVertex):
    """
    Represents a vertex at the end of a polymer chain.

    Do not instantiate this class directly, use :class:`.HeadVertex`
    or :class:`.TailVertex` instead.

    """

    def place_building_block(self, building_block, edges):
        if building_block.get_num_functional_groups() != 1:
            return super().place_building_block(building_block, edges)

        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg, = building_block.get_functional_groups()
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
                'The building block of a polymer '
                'must have 1 or 2 functional groups.'
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
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }
