"""
Metal Complex Vertices
======================

"""

from scipy.spatial.distance import euclidean

from stk.utilities import get_projection

from ..topology_graph import Vertex


class MetalVertex(Vertex):
    """
    Places the metal in a :class:`.MetalComplex`.

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


class MonoDentateLigandVertex(Vertex):
    """
    Places monodentate ligand in a :class:`.MetalComplex`.

    """

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        assert (
            building_block.get_num_functional_groups() == 1
        ), (
            f'{building_block} needs to have exactly 1 functional '
            'group but has '
            f'{building_block.get_num_functional_groups()}.'
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

        return {0: edges[0].get_id()}


class BiDentateLigandVertex(Vertex):
    """
    Places bidentate ligand in a :class:`.MetalComplex`.

    """

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        assert (
            building_block.get_num_functional_groups() == 2
        ), (
            f'{building_block} needs to have exactly 2 functional '
            'groups but has '
            f'{building_block.get_num_functional_groups()}.'
        )

        fg0_position, fg1_position = (
            building_block.get_centroid(fg.get_placer_ids())
            for fg in building_block.get_functional_groups()
        )
        edge_position1, edge_position2 = (
            edge.get_position() for edge in edges
        )
        building_block = building_block.with_rotation_between_vectors(
            start=fg1_position-fg0_position,
            target=edge_position2-edge_position1,
            origin=building_block.get_centroid(),
        )

        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        core_to_placer = placer_centroid - core_centroid

        fg0_position, fg1_position = (
            building_block.get_centroid(fg.get_placer_ids())
            for fg in building_block.get_functional_groups()
        )
        fg_vector = fg1_position - fg0_position

        fg_vector_projection = get_projection(
            start=core_to_placer,
            target=fg_vector,
        )

        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        building_block = building_block.with_rotation_between_vectors(
            start=core_to_placer - fg_vector_projection,
            target=edge_centroid - self._position,
            origin=building_block.get_centroid(),
        )

        return building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
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
