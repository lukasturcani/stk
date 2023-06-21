import typing

T = typing.TypeVar("T")


class MutationRecord(typing.Generic[T]):
    """
    Abstract base class for a record of a mutation operation.

    Notes:

        You might notice that the public methods of this abstract base
        class are implemented. This is just a default implementation, which
        can be used directly by users and subclasses, but can also be
        freely replaced during subclass implementation, if need be.

    """

    def __init__(self, molecule_record: T, mutator_name: str) -> None:
        """
        Parameters:
            molecule_record:
                The molecule produced by the mutation operation.

            mutator_name:
                The name of the mutator which carried out the mutator.

        """

        self._molecule_record = molecule_record
        self._mutator_name = mutator_name

    def get_molecule_record(self) -> T:
        """
        Get the molecule record produced by the mutation.

        Returns:
            The molecule record.

        """

        return self._molecule_record

    def get_mutator_name(self) -> str:
        """
        Get the name of the mutator which carried out the mutation.

        Returns:
            The name of the mutator.

        """

        return self._mutator_name
