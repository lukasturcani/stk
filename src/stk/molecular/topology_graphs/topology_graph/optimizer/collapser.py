"""
Collapser
=========

"""

from .optimizer import Optimizer
import mchammer as mch


class Collapser(Optimizer):
    """
    Perform rigid-body collapse of molecule.

    Implementation: [1]_

    .. [1] https://github.com/andrewtarzia/MCHammer.

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
        return state.with_position_matrix(
            mch_result.get_final_position_matrix()
        )

