class EdgeGroup:
    def __init__(self, edges):
        self._edge_ids = tuple(edge.get_id() for edge in edges)

    def get_edge_ids(self):
        yield from self._edge_ids
