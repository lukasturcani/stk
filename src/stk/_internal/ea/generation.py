import typing
from collections.abc import Iterable, Iterator
from dataclasses import dataclass

from stk._internal.ea.crossover.record import CrossoverRecord
from stk._internal.ea.mutation.record import MutationRecord


@dataclass(frozen=True, slots=True)
class FitnessValues:
    raw: typing.Any
    normalized: float


T = typing.TypeVar("T")


class Generation(typing.Generic[T]):
    """
    An abstract base class for EA generations.

    Notes:
        You might notice that the public methods of this abstract base
        class are implemented. This is just a default implementation, which
        you can ignore or override when implementing subclasses.
    """

    def __init__(
        self,
        fitness_values: dict[T, FitnessValues],
        mutation_records: Iterable[MutationRecord[T]],
        crossover_records: Iterable[CrossoverRecord[T]],
    ) -> None:
        """
        Parameters:
            fitness_values:
                The records of molecules in the generation.
            mutation_records (list[MutationRecord[T]]):
                The records of mutations done during the generation.
            crossover_record (list[CrossoverRecord[T]]):
                The records of crossover operations done during the
                generation.
        """
        self._fitness_values = dict(fitness_values)
        self._mutation_records = tuple(mutation_records)
        self._crossover_records = tuple(crossover_records)

    def get_molecule_records(self) -> Iterator[T]:
        """
        Yield the molecule records in the generation.

        Yields:
            A molecule record.
        """
        yield from self._fitness_values

    def get_mutation_records(self) -> Iterator[MutationRecord[T]]:
        """
        Yield the mutation records in the generation.

        Yields:
            A mutation record.
        """
        yield from self._mutation_records

    def get_crossover_records(self) -> Iterator[CrossoverRecord[T]]:
        """
        Yield the crossover records in the generation.

        Yields:
            A crossover record.
        """
        yield from self._crossover_records
