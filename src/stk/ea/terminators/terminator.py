class Terminator:
    """
    Checks if the exit criterion for the EA has been satisfied.

    """

    def terminate(self, progress):
        """
        Check to see the the EA should stop.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a
            subpopulation.

        Returns
        -------
        :class:`bool`
            ``True`` if the EA should stop and ``False`` otherwise.

        """

        raise NotImplementedError()
