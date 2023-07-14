import typing
from collections import Counter
from collections.abc import Iterable, Iterator

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.key_makers.molecule import MoleculeKeyMaker

BatchKey: typing.TypeAlias = frozenset[tuple[str, int]]

T = typing.TypeVar("T", bound=MoleculeRecord)


class Batch(typing.Generic[T]):
    """
    Represents a batch of molecule records.

    Batches can be compared, the comparison is based on their
    fitness values. Batches can also be iterated through, this
    iterates through all the records in the batch.

    Examples:

        *Sorting Batches by Fitness Value*

        Sorting batches causes them to be sorted by fitness value.

        .. testcode:: sorting-batches-by-fitness-value

            import stk

            record1 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            record2 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            record3 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )

            batches = (
                stk.Batch(
                    records={record1: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records={record2: 2},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records={record3: 3},
                    key_maker=stk.Inchi(),
                ),
            )
            sorted_batches = sorted(batches)

        .. testcode:: sorting-batches-by-fitness-value
            :hide:

            assert (
                sorted_batches[0] < sorted_batches[1] < sorted_batches[2]
            )

        *Comparing Batches by Fitness Value*

        Comparison is also based on fitness value

        .. testcode:: comparing-batches-by-fitness-value

            import stk

            record1 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            batch1 = stk.Batch(
                records={record1: 1},
                key_maker=stk.Inchi(),
            )

            record2 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            batch2 = stk.Batch(
                records={record2: 2},
                key_maker=stk.Inchi(),
            )

            if batch1 < batch2:
                print('batch1 has a smaller fitness value than batch2.')

        .. testoutput:: comparing-batches-by-fitness-value
            :hide:

            batch1 has a smaller fitness value than batch2.

        *Iterating Through Molecule Records in a Batch*

        Batches can be iterated through to get the molecule records in the
        batch

        .. testcode:: iterating-through-molecule-records-in-a-batch

            import stk

            record1 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            record2 = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=[
                        stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                    ],
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            batch = stk.Batch(
                records={record1: 1, record2: 2},
                key_maker=stk.Inchi(),
            )
            for record in batch:
                # Do stuff with record.
                pass
    """

    __slots__ = "_records", "_fitness_value", "_identity_key"

    def __init__(
        self,
        records: dict[T, float] | Iterable[tuple[T, float]],
        key_maker: MoleculeKeyMaker,
    ) -> None:
        """
        Parameters:
            records (dict[T, float]):
                The records which form the batch, mapped to their
                fitness values.
            key_maker:
                Used to make keys for molecules, which are used to
                determine the identity key of the batch. If two
                batches have the same molecule keys, the same number of
                times, they will have the same identity key.
        """
        self._records = dict(records)
        self._fitness_value = sum(self._records.values())
        molecules = (record.get_molecule() for record in self._records)
        self._identity_key = frozenset(
            Counter(map(key_maker.get_key, molecules)).items()
        )

    def get_size(self) -> int:
        """
        Get the number of molecules in the batch.

        Returns:
            The number of molecules in the batch.
        """
        return len(self._records)

    def get_fitness_value(self) -> float:
        """
        Get the fitness value of the batch.

        Returns:
            The fitness value.
        """
        return self._fitness_value

    def get_identity_key(self) -> BatchKey:
        """
        Get the identity key of the batch.

        If two batches hold the same molecules, the same number of
        times, they will have the same identity key.

        Returns:
            A hashable object which can be used to compare if two
            batches have the same identity.
        """
        return self._identity_key

    def __iter__(self) -> Iterator[T]:
        return iter(self._records)

    def __eq__(self, other: typing.Any) -> bool:
        if not isinstance(other, Batch):
            return NotImplemented
        return self._fitness_value == other._fitness_value

    def __gt__(self, other: "Batch[T]") -> bool:
        return self._fitness_value > other._fitness_value

    def __ge__(self, other: "Batch[T]") -> bool:
        return self._fitness_value >= other._fitness_value

    def __lt__(self, other: "Batch[T]") -> bool:
        return self._fitness_value < other._fitness_value

    def __le__(self, other: "Batch[T]") -> bool:
        return self._fitness_value <= other._fitness_value

    def __repr__(self) -> str:
        return f"Batch({self._fitness_value})"

    def __str__(self) -> str:
        return repr(self)
