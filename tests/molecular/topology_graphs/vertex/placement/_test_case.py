class _TestCase:
    def __init__(
        self,
        vertex,
        edges,
        building_block,
        position,
        alignment_tests,
        functional_group_edges,
    ):
        self.vertex = vertex
        self.edges = edges
        self.building_block = building_block
        self.position = position
        self.alignment_tests = alignment_tests
        self.functional_group_edges = functional_group_edges

    def __str__(self):
        return (
            f'TestCase({self.vertex}, {self.edges}, '
            f'{self.building_block}, {self.position}, '
            f'{self.alignment_tests}, {self.functional_group_edges})'
        )

    def __repr__(self):
        return str(self)
