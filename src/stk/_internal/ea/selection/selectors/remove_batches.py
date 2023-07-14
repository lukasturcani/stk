import typing
from collections.abc import Iterator

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import Batch

from .selector import ExcludedBatches, IncludedBatches, Selector

T = typing.TypeVar("T", bound=MoleculeRecord)


class RemoveBatches(Selector[T]):
    """
    Prevents a :class:`.Selector` from selecting some batches.

    Examples:
        *Removing Batches From Selection*

        You want to use :class:`.Roulette` selection on all but the
        5 :class:`.Worst` batches

        .. testcode:: removing-batches-from-selection

            import stk

            selector = stk.RemoveBatches(
                remover=stk.Worst(5),
                selector=stk.Roulette(20),
            )
            population = {
                stk.MoleculeRecord(
                    topology_graph=stk.polymer.Linear(
                        building_blocks=[
                            stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                        ],
                        repeating_unit='A',
                        num_repeating_units=2,
                    ),
                ): i
                for i in range(100)
            }
            for batch in selector.select(population):
                # Do stuff with batch. It was selected with roulette
                # selection and is not one of the worst 5 batches.
                pass
    """

    def __init__(self, remover: Selector[T], selector: Selector[T]) -> None:
        """
        Parameters:
            remover : :class:`.Selector`
                Selects batches of molecules, which cannot be yielded by
                `selector`.

            selector : :class:`.Selector`
                Selects batches of molecules, except those selected by
                `remover`.
        """
        self._remover = remover
        self._selector = selector

    def select(
        self,
        population: dict[T, float],
        included_batches: "IncludedBatches" = None,
        excluded_batches: "ExcludedBatches" = None,
    ) -> Iterator[Batch[T]]:
        removed_batches = {
            batch.get_identity_key()
            for batch in self._remover.select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        }
        if excluded_batches is not None:
            removed_batches |= excluded_batches

        yield from self._selector.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=removed_batches,
        )
