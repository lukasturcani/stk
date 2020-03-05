from ..topology_graph import Edge


class _CofEdge(Edge):
    def __init__(
        self,
        parent_id,
        id,
        vertex1,
        vertex2,
        periodicity=(0, 0, 0),
        position=None,
    ):
        super().__init__(id, vertex1, vertex2, periodicity, position)
        self._parent_id = parent_id

    def get_parent_id(self):
        return self._parent_id

    def clone(self):
        obj = super().clone()
        obj._parent_id = self._parent_id
        return obj
