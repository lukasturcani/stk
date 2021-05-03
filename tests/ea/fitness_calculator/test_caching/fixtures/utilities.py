class Counter:
    """
    A counter.

    """

    def __init__(self):
        self._count = 0

    def get_count(self, *args, **kwargs):
        """
        Return the number of times this method was called.

        Returns
        -------
        :class:`int`
            The number of times the method was called.

        """

        self._count += 1
        return self._count
