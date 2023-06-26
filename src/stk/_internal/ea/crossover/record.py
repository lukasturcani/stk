import typing

T = typing.TypeVar("T")


class CrossoverRecord(typing.Generic[T]):
    """
    Abstract base class for a record of a crossover operation.

    Notes:

        You might notice that the public methods of this abstract base
        class are implemented. This is just a default implementation, which
        can be used directly by users and subclasses, but can also be
        freely replaced during subclass implementation, if need be.

    """

    def __init__(self, molecule_record: T, crosser_name: str) -> None:
        """
        Parameters:

            molecule_record:
                The molecule produced by the crossover operation.

            crosser_name:
                The name of the crosser which carried out the crossover.

        """

        self._molecule_record = molecule_record
        self._crosser_name = crosser_name

    def get_molecule_record(self) -> T:
        """
        Get the molecule record produced by the crossover.

        Returns:
            The molecule record.

        """

        return self._molecule_record

    def get_crosser_name(self) -> str:
        """
        Get the name of the crosser which carried out the crossover.

        Returns:
            The name of the crosser.

        """

        return self._crosser_name
