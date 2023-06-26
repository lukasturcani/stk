import typing
from collections.abc import Iterator, Sequence

from stk._internal.ea.crossover.record import CrossoverRecord

T = typing.TypeVar("T")


class MoleculeCrosser(typing.Protocol[T]):
    """
    Performs crossover operations.

    Crossers take multiple molecules and recombine them to make
    new, offspring, molecules.

    Examples:

        *Implementing crossover operations*

        You only need to implement :meth:`.cross`. The source code of any
        of the classes listed in :mod:`.crosser` can serve as good
        examples.

    """

    def cross(self, records: Sequence[T]) -> Iterator[CrossoverRecord[T]]:
        """
        Cross `records`.

        Parameters:
            records (list[T]):
                The molecule records on which a crossover operation is
                performed.

        Yields:
            A record of a crossover operation.

        """
        ...
