
class AnyTerminator(Terminator):
    """
    Checks if any :class:`Terminator` has satisfied its exit condition.

    """

    def __init__(self, *Terminators):
        """
        Initialize a :class:`AnyTerminator` instance.

        Parameters
        ----------
        *Terminators : :class:`Terminator`
            :class:`Terminator` objects which are checked to see if their
            exit conditions have been satisfied.

        """

        self._Terminators = Terminators

    def terminate(self, progress):
        """
        Check to see if any exit condition has been satisfied.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if any :class:`Terminator` in :attr:`Terminators` has
            satisfied its exit condition.

        """

        return any(
            Terminator.terminate(progress)
            for Terminator in self._Terminators
        )


