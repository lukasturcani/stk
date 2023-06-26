import typing

from stk._internal.ea.mutation.record import MutationRecord

T = typing.TypeVar("T")


class MoleculeMutator(typing.Protocol[T]):
    """
    Performs mutation operations.

    Examples:

        *Implementation mutation operations*

        You only need to implement :meth:`.mutate`. The source code of any
        of the classes listed in :mod:`.mutator` can serve as good
        examples.

    """

    def mutate(self, record: T) -> MutationRecord[T] | None:
        """
        Return a mutant of `record`.

        Parameters:
            record:
                The molecule to be mutated.

        Returns:
            A record of the mutation or ``None``
            if `record` cannot be mutated.
        """
        ...
