"""
MCHammer
========

"""


import mchammer as mch


class MCHammer(Optimizer):
    """
    Perform Monte Carlo based optimisation of long-bonds in molecule.

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

