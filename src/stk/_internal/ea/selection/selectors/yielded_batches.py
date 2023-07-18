import typing

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import Batch, BatchKey
from stk._internal.key_makers.molecule import MoleculeKeyMaker

T = typing.TypeVar("T", bound=MoleculeRecord)


class YieldedBatches(typing.Generic[T]):
    """
    Keeps track of batches yielded by :meth:`.Selector.select`.
    """

    __slots__ = ("_molecules", "_batches", "_num", "_key_maker")

    def __init__(self, key_maker: MoleculeKeyMaker) -> None:
        self._key_maker = key_maker

        # Has all molecules yielded by select().
        self._molecules: set[str] = set()
        # Has the identity_key() of all batches yielded by select().
        self._batches: set[BatchKey] = set()
        # Counts the total number of times select() has yielded.
        self._num = 0

    def update(self, batch: Batch[T]) -> typing.Self:
        """
        Update tracked data with a new `batch`.

        Parameters:
            batch:
                A batch yielded by :meth:`.Selector.select`.
        Returns:
            The data tracker.
        """
        molecules = (record.get_molecule() for record in batch)
        self._molecules.update(map(self._key_maker.get_key, molecules))
        self._batches.add(batch.get_identity_key())
        self._num += 1
        return self

    def get_num(self) -> int:
        """
        Get the number of times :meth:`.Selector.select` has yielded.

        Returns:
            The total number of times :meth:`.Selector.select` has
            yielded.
        """
        return self._num

    def is_yielded_batch(self, batch: Batch[T]) -> bool:
        """
        Check if `batch` has already been yielded.

        Parameters:
            batch:
                The batch to check.
        Returns:
            ``True`` if `batch` has already been yielded.
        """
        return batch.get_identity_key() in self._batches

    def is_unyielded_batch(self, batch: Batch[T]) -> bool:
        """
        Check if `batch` has not been yielded.

        Parameters:
            batch:
                The batch to check.
        Returns:
            ``True`` if `batch` has not been yielded.
        """
        return batch.get_identity_key() not in self._batches

    def has_yielded_molecules(self, batch: Batch[T]) -> bool:
        """
        Check if `batch` contains any previously yielded molecules.

        Parameters:
            batch:
                The batch to check.
        Returns:
            ``True`` if `batch` contains any molecules which have
            previously been yielded.
        """
        return any(
            self._key_maker.get_key(molecule) in self._molecules
            for molecule in (record.get_molecule() for record in batch)
        )

    def has_no_yielded_molecules(self, batch: Batch[T]) -> bool:
        """
        Check if `batch` consists only of unyielded molecules.

        Parameters:
            batch:
                The batch to check.
        Returns:
            ``True`` if `batch` does not have any previously yielded
            molecules.
        """
        return all(
            self._key_maker.get_key(molecule) not in self._molecules
            for molecule in (record.get_molecule() for record in batch)
        )
