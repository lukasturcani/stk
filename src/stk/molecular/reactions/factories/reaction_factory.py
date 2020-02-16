class ReactionFactory:
    """
    Creates :class:`.Reaction` instances.

    """

    def get_reaction(self, functional_groups):
        """
        Get a reaction to use on the `functional_groups`.

        Parameters
        ----------
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
