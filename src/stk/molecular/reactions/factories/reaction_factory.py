class ReactionFactory:
    """
    Creates :class:`.Reaction` instances.

    """

    def get_reaction(
        self,
        construction_state,
        edges,
        edge_group,
        functional_groups,
    ):
        """
        Get a reaction to use on the `functional_groups`.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The state of the current construction.

        edges : :class:`iterable` of :class:`.Edge`
            The edges on which the `functional_groups` are found.

        edge_group : :class:`.EdgeGroup`
            The edge group on which the `functional_groups` are found.

        functional_groups : :class:`iterable`
            An :class:`iterable` holding :class:`.FunctionalGroup`
            instances, for which a :class:`.Reaction` should be
            returned.

        Returns
        -------
        :class:`.Reaction`
            The reaction to use on the `functional_groups`.

        """

        raise NotImplementedError()
