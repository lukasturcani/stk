class ReactionFactory:
    """
    Creates :class:`.Reaction` instances.

    """

    def get_reaction(self, construction_state, edge_group):
        """
        Get a reaction to use on the `functional_groups`.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The state of the current construction.

        edge_group : :class:`.EdgeGroup`
            The edge group for which a reaction should be found.

        Returns
        -------
        :class:`.Reaction`
            The reaction to use on the `edge_group`.

        """

        raise NotImplementedError()
