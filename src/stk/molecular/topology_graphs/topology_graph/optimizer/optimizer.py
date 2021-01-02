class Optimizer:
    """
    An abstract base class for optimizers.

    An optimizer is used to change the structure of molecule under
    construction to be more realistic.

    """

    def optimize(self, state):
        """
        Optimize the structure of a molecule under construction.

        Parameters
        ----------
        state : :class:`.ConstructionState`
            The molecule being constructed.

        Returns
        -------
        :class:`.ConstructionState`
            The optimized construction state.

        """

        raise NotImplementedError()
