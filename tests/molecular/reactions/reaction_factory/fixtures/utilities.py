import stk


class MockConstructionState(stk.ConstructionState):
    def __init__(self, position_matrix=None):
        self.position_matrix = position_matrix

    def get_position_matrix(self):
        return self.position_matrix


class MockEdge(stk.Edge):
    def __init__(self, periodicity):
        self.periodicity = periodicity

    def get_periodicity(self):
        return self.periodicity
