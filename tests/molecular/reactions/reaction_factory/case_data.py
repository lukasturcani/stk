class CaseData:
    """
    A :class:`.ReactionFactory` test case.

    Attributes
    ----------
    factory : :class:`.ReactionFactory`
        The reaction factory to test.

    construction_state : :class:`.ConstructionState`
        The construction state to pass to the factory.

    edge_group : :class:`.EdgeGroup`
        The edge group to pass to the factory.

    reaction_result : :class:`.ReactionResult`
        The expected result of the reaction returned by the factory.

    """

    def __init__(
        self,
        factory,
        construction_state,
        edge_group,
        reaction_result,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        factory : :class:`.ReactionFactory`
            The reaction factory to test.

        construction_state : :class:`.ConstructionState`
            The construction state to pass to the factory.

        edge_group : :class:`.EdgeGroup`
            The edge group to pass to the factory.

        reaction_result : :class:`.ReactionResult`
            The expected result of the reaction returned by the
            factory.

        """

        self.factory = factory
        self.construction_state = construction_state
        self.edge_group = edge_group
        self.reaction_result = reaction_result
