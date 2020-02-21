class _TestCase:
    def __init__(
        self,
        factory,
        construction_state,
        edges,
        edge_group,
        functional_groups,
        reaction_result,
    ):
        self.factory = factory
        self.construction_state = construction_state
        self.edges = edges
        self.edge_group = edge_group
        self.functional_groups = functional_groups
        self.reaction_result = reaction_result
