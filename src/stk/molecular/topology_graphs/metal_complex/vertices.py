import numpy as np
from scipy.spatial.distance import euclidean

from ..topology_graph import Vertex


class _MetalVertex(Vertex):
    """
    Places the metal in a :class:`.MetalComplex`.

    """

    def place_building_block(self, building_block, edges):
        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        edges = sorted(edges, key=lambda i: i.get_id())
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class _MonoDentateLigandVertex(Vertex):
    """
    Places monodentate ligand in a :class:`.MetalComplex`.

    """

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg, = building_block.get_functional_groups(0)
        fg_centroid = building_block.get_centroid(
            atom_ids=fg.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        return building_block.with_rotation_between_vectors(
            start=fg_centroid - core_centroid,
            target=edge_centroid - self._position,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        edges = sorted(edges, key=lambda i: i.get_id())
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class _BiDentateLigandVertex(Vertex):
    """
    Places bidentate ligand in a :class:`.MetalComplex`.

    """

    def __init__(
        self,
        id,
        position,
    ):
        """
        Initialize a :class:`._BiDentateLigandVertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        """

        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        return clone

    def place_building_block(self, building_block, edges):
        # Translate building block to vertex position.
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )

        # Align vector between 2 edges with vector between centroid of
        # placers in 2 FGs.
        fg0, fg1 = building_block.get_functional_groups()
        fg0_position = building_block.get_centroid(
            atom_ids=fg0.get_placer_ids(),
        )
        fg1_position = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        start = fg1_position - fg0_position
        # Vector between connected edges.
        c_edge_positions = [
            i.get_position() for i in edges
        ]
        target = c_edge_positions[1] - c_edge_positions[0]
        building_block = building_block.with_rotation_between_vectors(
            start=start,
            target=target,
            origin=building_block.get_centroid(),
        )

        # Align vector between edge-self.position with vector between
        # placer centroid and core of the molecule centroid.
        # Importantly, we use a projection of the placer-core vector
        # that is orthogonal to the FG-FG vector in a bidentate
        # ligand for this alignment.
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )

        fg0, fg1 = building_block.get_functional_groups()
        fg0_position = building_block.get_centroid(
            atom_ids=fg0.get_placer_ids(),
        )
        fg1_position = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        fg_vector = fg1_position - fg0_position
        placer_to_core_vector = placer_centroid - core_centroid
        proj_onto_fg_vector = fg_vector * np.dot(
            placer_to_core_vector,
            fg_vector
        ) / np.dot(fg_vector, fg_vector)
        orthogonal_vector = proj_onto_fg_vector - placer_to_core_vector
        start = orthogonal_vector
        target = self._position - edge_centroid
        building_block = building_block.with_rotation_between_vectors(
            start=start,
            target=target,
            origin=building_block.get_centroid(),
        )

        # Translate building block to vertex position.
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }
