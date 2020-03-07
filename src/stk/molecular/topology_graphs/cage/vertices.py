import numpy as np
from scipy.spatial.distance import euclidean

from stk.utilities import (
    get_acute_vector,
    get_plane_normal,
    vector_angle,
    normalize_vector,
)
from ..topology_graph import Vertex


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
        self._use_neighbor_placement = use_neighbor_placement
        # The edge which is used to align the :class:`.BuildingBlock`
        # placed on the vertex. The first :class:`.FunctionalGroup`
        # is rotated such that it lies exactly on this :class:`.Edge`.
        # Must be between ``0`` and the number of edges the vertex is
        # connected to.
        self._aligner_edge = aligner_edge
        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        clone._use_neighbor_placement = self._use_neighbor_placement
        return clone

    def _with_aligner_edge(self, aligner_edge):
        self._aligner_edge = aligner_edge
        return self

    def with_aligner_edge(self, aligner_edge):
        return self.clone()._with_aligner_edge(aligner_edge)

    def use_neighbor_placement(self):
        return self._use_neighbor_placement

    @classmethod
    def init_at_center(cls, id, vertices):
        return cls(
            id=id,
            position=(
                sum(vertex.get_position() for vertex in vertices)
                / len(vertices)
            ),
        )

    def __str__(self):
        return (
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class _LinearCageVertex(_CageVertex):
    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        start = fg_centroid - self._position
        edge_coord = edges[self._aligner_edge].get_position()
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        building_block = building_block.with_rotation_between_vectors(
            start=start,
            target=edge_coord - edge_centroid,
            origin=self._position,
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=core_centroid - self._position,
                target=self._position,
                axis=normalize_vector(
                    edges[0].get_position() - edges[1].get_position()
                ),
                origin=self._position,
            )
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg = next(building_block.get_functional_groups(0))
        fg0_position = building_block.get_centroid(fg.get_placer_ids())
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(sorted(
                edges,
                key=lambda edge: euclidean(
                    edge.get_position(),
                    fg0_position,
                ),
            ))
        }


class _NonLinearCageVertex(_CageVertex):
    def place_building_block(self, building_block, edges):
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
        start = fg_bonder_centroid - self._position
        edge_coord = edges[self._aligner_edge].get_position()
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=start,
                target=edge_coord - edge_centroid,
                axis=edge_normal,
                origin=self._position,
            )
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        # The idea is to order the functional groups in building_block
        # by their angle with the vector running from the placer
        # centroid to fg0, going in the clockwise direction.
        # The edges are also ordered by their angle with the vector
        # running from the edge centroid to the aligner_edge,
        # going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg0_position = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
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
        functional_groups = sorted(
            range(building_block.get_num_functional_groups()),
            key=self._get_functional_group_angle(
                building_block=building_block,
                fg0_direction=fg0_direction,
                centroid=placer_centroid,
                axis=axis,
            ),
        )
        edges = sorted(edges, key=self._get_edge_angle(axis, edges))
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(functional_groups, edges)
        }

    def _get_functional_group_angle(
        self,
        building_block,
        fg0_direction,
        centroid,
        axis,
    ):

        def angle(fg_id):
            fg = next(building_block.get_functional_groups(fg_id))
            position = building_block.get_centroid(
                atom_ids=fg.get_placer_ids(),
            )
            fg_direction = position - centroid
            theta = vector_angle(fg0_direction, fg_direction)

            projection = fg_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def _get_edge_angle(self, axis, edges):
        aligner_edge = edges[self._aligner_edge]
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        aligner_edge_direction = (
            aligner_edge.get_position() - edge_centroid
        )

        def angle(edge):
            edge_direction = edge.get_position() - edge_centroid
            theta = vector_angle(
                vector1=edge_direction,
                vector2=aligner_edge_direction,
            )

            projection = edge_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle
