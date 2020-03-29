from .terminator import Terminator


class NumGenerations(Terminator):
    """
    Stop the EA after a certain number of generations.

    """

    def __init__(self, num_generations):
        """
        Initialize a :class:`NumGenerations` instance.

        Parameters
        ----------
        num_generations : :class:`int`
            The number of generations after which the EA should stop.

        """

        self._num_generations = num_generations

    def terminate(self, progress):
        """
        Check if a number of generations has passed.

        Parameters
        ----------
        progress : :class:`.Population`
            A population where every generation is a subpopulation.

        Returns
        -------
        ``True`` if :attr:`num_generations` or more has passed.

        """

        return len(progress.subpopulations) >= self._num_generations
