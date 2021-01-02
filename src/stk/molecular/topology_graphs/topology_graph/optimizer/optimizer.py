from ..construction_state import ConstructionState


class Optimizer:
    """
    An abstract base class for optimizers.

    An optimizer is used to change the structure of molecule under
    construction to be more realistic.

    """

    def optimize(self, state: ConstructionState) -> ConstructionState:
        """
        Optimize the structure of a molecule under construction.

        Parameters
        ----------
        state
            The molecule being constructed.

        Returns
        -------
        The optimized construction state.


        """

        raise NotImplementedError()
