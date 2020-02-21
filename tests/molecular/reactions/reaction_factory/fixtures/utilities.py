import stk


class MockConstructionState(stk.ConstructionState):
    def __init__(
        self,
        edges,
        edge_functional_groups,
        position_matrix=None,
    ):
        self._edges = edges
        self._edge_functional_groups = edge_functional_groups
        self._position_matrix = position_matrix


class MockEdge(stk.Edge):
    def __init__(self, id, periodicity):
        self._id = id
        self._periodicity = periodicity
