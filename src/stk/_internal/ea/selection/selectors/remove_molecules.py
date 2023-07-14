import typing
from collections.abc import Iterator

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import Batch

from .selector import ExcludedBatches, IncludedBatches, Selector

T = typing.TypeVar("T", bound=MoleculeRecord)


class RemoveMolecules(Selector[T]):
    """
    Prevents a :class:`.Selector` from selecting some molecules.

    Examples:

        *Removing Molecules From Selection*

        .. testcode:: removing-molecules-from-selection

            import stk

            selector = stk.RemoveMolecules(
                # Do not select any molecules from the top 5 batches
                # of size 3.
                remover=stk.Best(num_batches=5, batch_size=3),
                # Select the ret of the molecules using roulette.
                selector=stk.Roulette(num_batches=20, batch_size=3),
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
                for i in range(10)
            }
            for batch in selector.select(population):
                # Do stuff with batch. The batch is guaranteed not to
                # contain any molecules which are found in the best 5
                # batches of size 3.
                pass
    """

    def __init__(self, remover: Selector[T], selector: Selector[T]) -> None:
        """
        Parameters:
            remover:
                Selects batches molecules, any molecule selected cannot
                be selected by `selector`.

            selector:
                Selects batches of molecules, not containing any molecules
                selected by `remover`.
        """
        self._remover = remover
        self._selector = selector

    def select(
        self,
        population: dict[T, float],
        included_batches: "IncludedBatches" = None,
        excluded_batches: "ExcludedBatches" = None,
    ) -> Iterator[Batch[T]]:
        remover_batches = self._remover.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        removed = {record for batch in remover_batches for record in batch}
        valid_population = {
            record: population[record]
            for record in population
            if record not in removed
        }
        yield from self._selector.select(
            population=valid_population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
