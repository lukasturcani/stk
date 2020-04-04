"""
All Terminators
===============

"""

from .terminator import Terminator


class AllTerminators(Terminator):
    """
    Makes sure multiple terminators agree to terminate.

    This is a compound :class:`.Terminator` it agree to terminate an
    EA only if all of the its component :class:`.Terminator` instances
    agree to do so.

    Examples
    --------
    *Making Sure Multiple Termination Conditions Have Been Reached*

    Sometimes you want to make sure multiple termination conditions
    have been reached. For example, that a specific number of
    generation have passed, and that the fitness has plateaued.

    .. code-block:: python

        import stk

        terminator = stk.AllTerminators(
            # Only terminate if the number of generations is bigger
            # than or equal to 50 and the fittest member of the
            # population have not changed for the last 3 generations.
            terminators=(
                stk.NumGenerations(50),
                stk.FitnessPlateau(3),
            ),
        )
        # Use as normal.
        should_terminate = terminator.terminate(population)

    """

    def __init__(self, terminators):
        """
        Initialize an :class:`.AllTerminators` instance.

        Parameters
        ----------
        terminators : :class:`tuple` of :class:`.Terminator`
            :class:`.Terminator` instances, each of  which is checked
            to see if its exit condition have been satisfied.

        """

        self._terminators = terminators

    def terminate(self, population):
        return all(
            terminator.terminate(population)
            for terminator in self._terminators
        )
