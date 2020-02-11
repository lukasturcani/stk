class _ConstructionState:
    """

    """

    def __init__(self, vertex_assignments, vertex_edges):
        self._vertex_assignments = {
            v.get_id(): bb for v, bb in vertex_assignments.items()
        }
        self._vertices = {
            v.get_id(): v for v in vertex_assignments
        }
        self._vertex_edges = dict(vertex_edges)

    def with_placement_results(self, results):
        pass

    def get_vertex_assignments(self):
        """

        """

        return dict(self._vertex_assignments)

    def get_vertices(self, vertex_ids):
        """

        """

        for vertex_id in vertex_ids:
            yield self._vertices[vertex_id]

    def get_vertex_edges(self, vertex_ids):
        pass

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
        _ConstructionState.__init__(
            self=clone,
            vertex_assignments={
                vertex: self._vertex_assignments[id_]
                for id_, vertex in self._vertices.items()
            },
        )
        return clone
