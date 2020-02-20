class ConstructionState:
    """

    """

    def __init__(self, building_block_vertices, vertex_edges):
        """

        """

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
        self._vertex_edges = dict(vertex_edges)

    def with_placement_results(self, results):
        pass

    def get_vertex_building_block(self, vertex_id):
        """

        """

        return self._vertex_building_blocks[vertex_id]

    def get_vertices(self, vertex_ids):
        """

        """

        for vertex_id in vertex_ids:
            yield self._vertices[vertex_id]

    def get_vertex_edges(self, vertex_ids):
        for vertex_id in vertex_ids:
            yield self._vertex_edges[vertex_id]

    def get_edge_functional_groups(self):
        # Returns a dictionary.
        pass

    def with_reaction_results(self, results):
        pass

    def _with_vertices(self, vertices):
        self._vertices = {
            vertex.get_id(): vertex for vertex in vertices
        }
        return self

    def with_vertices(self, vertices):
        return self.clone()._with_vertices(vertices)

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        ConstructionState.__init__(
            self=clone,
            vertex_assignments={
                vertex: self._vertex_assignments[id_]
                for id_, vertex in self._vertices.items()
            },
        )
        return clone

    def get_position_matrix(self):
        """

        """

        pass

    def get_atoms(self):
        pass

    def get_bonds(self):
        pass

    def get_atom_infos(self):
        pass

    def get_reaction_infos(self):
        pass

    def get_building_block_counts(self):
        pass
