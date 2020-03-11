import numpy as np
from collections import defaultdict


class _GraphState:
    def __init__(
        self,
        building_block_vertices,
        edges,
        lattice_constants,
    ):
        self._vertex_building_blocks = {
            vertex.get_id(): building_block
            for building_block, vertices
            in building_block_vertices.items()
            for vertex in vertices
        }
        self._vertices = {
                vertex.get_id(): vertex
                for vertices in building_block_vertices.values()
                for vertex in vertices
        }
        self._edges = edges
        self._lattice_constants = lattice_constants
        self._vertex_edges = self._get_vertex_edges()

    def _get_vertex_edges(self):
        vertex_edges = defaultdict(list)
        for edge in self._edges:
            vertex_ids = (edge.get_vertex1_id(), edge.get_vertex2_id())

            if edge.is_periodic():
                for vertex_id in vertex_ids:
                    periodic_edge = self._get_periodic_edge(
                        edge=edge,
                        reference=vertex_id,
                    )
                    vertex_edges[vertex_id].append(periodic_edge)
            else:
                for vertex_id in vertex_ids:
                    vertex_edges[vertex_id].append(edge)
        return vertex_edges

    def _get_periodic_edge(self, edge, reference):
        vertex1 = self._vertices[reference]
        id1, id2 = edge.get_vertex_ids()
        vertex2 = self._vertices[id1 if reference == id2 else id2]

        direction = 1 if reference == id1 else -1
        periodicity = np.array(edge.get_periodicity())
        end_cell = vertex1.get_cell() + direction*periodicity
        cell_shift = end_cell - vertex2.get_cell()

        shift = sum(
            axis_shift*constant
            for axis_shift, constant in zip(
                cell_shift,
                self._lattice_constants,
            )
        )
        position = (
            (vertex2.get_position()+shift+vertex1.get_position()) / 2
        )
        return edge.with_position(position)

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        clone._vertex_building_blocks = dict(
            self._vertex_building_blocks
        )
        clone._vertices = dict(self._vertices)
        clone._vertex_edges = dict(self._vertex_edges)
        clone._edges = self._edges
        clone._position_matrix = np.array(self._position_matrix)
        clone._atoms = list(self._atoms)
        clone._atom_infos = list(self._atom_infos)
        clone._bonds = list(self._bonds)
        clone._bond_infos = list(self._bond_infos)
        clone._edge_functional_groups = defaultdict.copy(
            self._edge_functional_groups
        )
        return clone

    def get_building_block(self, vertex_id):
        """

        """

        return self._vertex_building_blocks[vertex_id]

    def get_vertex(self, vertex_id):
        """

        """

        return self._vertices[vertex_id]

    def get_num_vertices(self):
        return len(self._vertices)

    def get_edge(self, edge_id):
        return self._edges[edge_id]

    def get_num_edges(self):
        return len(self._edges)

    def get_edges(self, vertex_id):
        return self._vertex_edges[vertex_id]

    def _with_vertices(self, vertices):
        self._vertices = {
            vertex.get_id(): vertex for vertex in vertices
        }
        return self

    def with_vertices(self, vertices):
        clone = self.clone()
        clone._graph_state = self._graph_state.with_vertices(vertices)
        return self.clone()._with_vertices(vertices)
