import numpy as np
from functools import partial
from scipy.spatial.distance import euclidean
from stk.utilities import vector_angle

from ...topology_graph import Vertex


class _CofVertex(Vertex):
    def __init__(self, id, position, aligner_edge):
        super().__init__(id, position)
        self._aligner_edge = aligner_edge

    @classmethod
    def init_at_center(cls, id, vertices, aligner_edge):
        obj = cls.__new__(cls)
        obj._id = id
        obj._position = np.array([0, 0, 0], dtype=np.float64)
        for count, vertex in enumerate(vertices, 1):
            obj._position += vertex.get_position()
        obj._position /= count
        obj._aligner_edge = aligner_edge
        return obj

    @classmethod
    def init_at_shifted_center(
        cls,
        vertices,
        cell_shifts,
        lattice_constants,
        aligner_edge,
    ):
        positions = []
        for vertex, cell_shift in zip(vertices, cell_shifts):
            shift = sum(
                dim_shift*constant
                for dim_shift, constant
                in zip(cell_shift, lattice_constants)
            )
            positions.append(vertex.get_position() + shift)

        position = np.divide(
            np.sum(positions, axis=0),
            len(positions),
        )
        return cls(id, position, aligner_edge)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        return clone

    def __str__(self):
        return (
            'Vertex(id={self._id}, '
            'position={self._position.tolist()}, '
            'aligner_edge={self._aligner_edge})'
        )


class _LinearCofVertex(_CofVertex):
    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg = next(building_block.get_functional_groups(0))
        fg_centroid = building_block.get_centroid(fg.get_placer_ids())
        start = fg_centroid - self._position
        target = edges[0].get_position() - edges[1].get_position()
        target *= 1 if self._aigner_edge == 0 else -1

        building_block = building_block.with_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._position,
        )
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=building_block.get_centroid() - self._position,
                target=[0, 0, 1],
                axis=target,
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
                key=partial(euclidean, fg0_position),
            ))
        }


class _NonLinearCofVertex(_CofVertex):
    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = _get_acute_vector(
            reference=self._position - building_block.get_centroid(),
            vector=normal,
        )
        building_block = building_block.with_rotation_between_vectors(
            start=normal,
            target=[0, 0, 1],
            origin=self._position,
        )
        fg = next(building_block.get_functional_groups(0))
        fg_centroid = building_block.get_centroid(fg.get_placer_ids())
        edge_coord = edges[self._aligner_edge].get_position()
        building_block = (
            building_block.with_rotation_to_minimize_angle(
                start=fg_centroid - self._position,
                target=edge_coord - self._position,
                axis=[0, 0, 1],
                origin=self._position,
            )
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg0 = next(building_block.get_functional_groups(0))
        fg0_position = building_block.get_centroid(
            atom_ids=fg0.get_placer_ids(),
        )
        centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        fg0_direction = fg0_position - centroid

        normal = building_block.get_plane_normal(
            atom_ids=building_block.get_placer_ids(),
        )
        normal = _get_acute_vector(
            reference=self._position - building_block.get_centroid(),
            vector=normal,
        )

        axis = np.cross(fg0_direction, normal)
        functional_groups = sorted(
            range(building_block.get_num_functional_groups()),
            key=self._get_functional_group_angle(
                building_block=building_block,
                fg0_direction=fg0_direction,
                centroid=centroid,
                axis=axis,
            ),
        )
        edges = sorted(edges, key=self._get_edge_angle(axis, edges))
        return {
            edge.get_id(): fg_id
            for edge, fg_id in zip(functional_groups, edges)
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
                return 2*np.pu - theta
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


def _get_acute_vector(reference, vector):
    if (
        # vector_angle is NaN if reference is [0, 0, 0].
        not np.allclose(reference, [0, 0, 0], atol=1e-5)
        and vector_angle(vector, reference) > np.pi/2
    ):
        return vector * -1
    return vector
