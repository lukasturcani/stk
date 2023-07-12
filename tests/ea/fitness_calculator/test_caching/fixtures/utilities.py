import typing


class Counter:
    def __init__(self) -> None:
        self._count = 0

    def get_count(self, *args: typing.Any, **kwargs: typing.Any) -> int:
        """
        Return the number of times this method was called.

        Returns:
            The number of times the method was called.
        """
        self._count += 1
        return self._count
