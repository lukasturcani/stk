from .implementation import Implementation


class Serial(Implementation):
    """
    A serial implementation of the default evolutionary algorithm.

    """

    def get_generations(self):
        yield from self._get_generations(map)
