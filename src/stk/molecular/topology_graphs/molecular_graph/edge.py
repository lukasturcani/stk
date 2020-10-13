class Edge:
    def __init__(self, functional_group1, functional_group2):
        self._functional_group1 = functional_group1
        self._functional_group2 = functional_group2

    def get_functional_group1(self):
        return self._functional_group1

    def get_functional_group2(self):
        return self._functional_group2
