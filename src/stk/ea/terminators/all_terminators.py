from .terminator import Terminator


class AllTerminators(Terminator):
    """
    Checks if all :class:`Terminator` objects return ``True``.

    """

    def __init__(self, *Terminators):
        """
        Initialize a :class:`AllTerminator` instance.

        Parameters
        ----------
        *Terminators : :class:`Terminator`
            :class:`Terminator` objects which are checked to see if their
            exit conditions have been satisfied.

        """

        self._Terminators = Terminators

    def terminate(self, progress):
        """
        Checks to see if all exit conditions have been satisfied.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if all :class:`Terminator` objects in :attr:`Terminators`
            have satisfied its exit condition.

        """

        return all(
            Terminator.terminate(progress)
            for Terminator in self._Terminators
        )
