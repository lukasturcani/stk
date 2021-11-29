"""
Cage Vertices
=============

"""

import numpy as np
from scipy.spatial.distance import euclidean

from stk.utilities import (
    get_acute_vector,
    get_plane_normal,
    normalize_vector,
)

from ..topology_graph import Vertex
from ..utilities import _EdgeSorter, _FunctionalGroupSorter


class _CageVertex(Vertex):
    """
    Represents a vertex of a :class:`.Cage`.

    """

    def __init__(
        self,
        id,
        position,
        use_neighbor_placement=True,
        aligner_edge=0,
    ):
        """
        Initialize a :class:`._CageVertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        use_neighbor_placement : :class:`bool`, optional
            If ``True``, the position of the vertex will be updated
            based on the neighboring functional groups.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        """

        self._use_neighbor_placement = use_neighbor_placement
        self._aligner_edge = aligner_edge
        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        clone._use_neighbor_placement = self._use_neighbor_placement
        return clone

    def _with_aligner_edge(self, aligner_edge):
        """
        Modify the instance.

        """

        self._aligner_edge = aligner_edge
        return self

    def with_aligner_edge(self, aligner_edge):
        """
        Return a clone with a different `aligner_edge`.

        Parameters
        ----------
        aligner_edge : :class:`int`
            The aligner edge of the clone.

        Returns
        -------
        :class:`._CageVertex`
            The clone. Has the same type as the original instance.

        """

        return self.clone()._with_aligner_edge(aligner_edge)

    def use_neighbor_placement(self):
        """
        ``True`` if the position should be updated based on neighbors.

        Returns
        -------
        :class:`bool`
            ``True`` if the position of the vertex should be updated
            based on the positions of functional groups on neighboring
            vertices.

        """

        return self._use_neighbor_placement

    @classmethod
    def init_at_center(cls, id, vertices):
        """
        Initialize a :class:`._CageVertex` in the middle of `vertices`.

        Parameters
        ----------
        id : :class:`int`
            The id of the initialized vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices at whose center this one needs to be.

        Returns
        -------
        :class:`._CageVertex`
            The new vertex.

        """

        return cls(
            id=id,
            position=(
                sum(vertex.get_position() for vertex in vertices)
                / len(vertices)
            ),
        )

    def get_aligner_edge(self):
        """
        Return the aligner edge of the vertex.

        Returns
        -------
        :class:`int`
            The aligner edge.

        """

        return self._aligner_edge

    def __str__(self):
        return (
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class LinearVertex(_CageVertex):
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
        fg_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        building_block = building_block.with_rotation_between_vectors(
            start=fg_centroid - self._position,
            target=edge_position - edge_centroid,
            origin=self._position,
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        return building_block.with_rotation_to_minimize_angle(
            start=core_centroid - self._position,
            target=self._position,
            axis=normalize_vector(
                edges[0].get_position() - edges[1].get_position()
            ),
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class NonLinearVertex(_CageVertex):
    def place_building_block(self, building_block, edges):
        assert (
            building_block.get_num_functional_groups() > 2
        ), (
            f'{building_block} needs to have more than 2 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        edge_normal = get_acute_vector(
            reference=edge_centroid,
            vector=get_plane_normal(
                points=np.array([
                    edge.get_position() for edge in edges
                ]),
            ),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        building_block = building_block.with_rotation_between_vectors(
            start=get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(
                    atom_ids=building_block.get_placer_ids(),
                ),
            ),
            target=edge_normal,
            origin=self._position,
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        return building_block.with_rotation_to_minimize_angle(
            start=fg_bonder_centroid - self._position,
            target=edge_position - edge_centroid,
            axis=edge_normal,
            origin=self._position,
        ).get_position_matrix()

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
            aligner_edge=edges[self._aligner_edge],
            axis=fg_sorter.get_axis(),
        )
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(
                fg_sorter.get_items(),
                edge_sorter.get_items(),
            )
        }


class UnaligningVertex(_CageVertex):
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

    @classmethod
    def init_at_center(cls, id, vertices):
        vertex = cls.__new__(cls)
        vertex._id = id
        vertex._position = (
            sum(vertex.get_position() for vertex in vertices)
            / len(vertices)
        )
        vertex._use_neighbor_placement = True
        vertex._aligner_edge = 0
        return vertex


class AngledVertex(_CageVertex):

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

        fg_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        building_block = building_block.with_rotation_between_vectors(
            start=fg_centroid - self._position,
            target=edge_position - edge_centroid,
            origin=self._position,
        )

        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        core_to_placer = placer_centroid - core_centroid
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )

        return building_block.with_rotation_between_vectors(
            start=core_to_placer,
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }
