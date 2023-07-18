import typing
from collections.abc import Iterator

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import Batch

from .selector import ExcludedBatches, IncludedBatches, Selector

T = typing.TypeVar("T", bound=MoleculeRecord)


class FilterMolecules(Selector[T]):
    """
    Allows a :class:`.Selector` to select only some molecules.

    Examples:

        *Using a Selection Algorithm on a Subset of Batches*

        .. testcode:: using-a-selection-algorithm-on-a-subset-of-batches

            import stk

            selector = stk.FilterMolecules(
                # Use only molecules in the top 5 batches of size 3.
                filter=stk.Best(num_batches=5, batch_size=3),
                # Select from those molecules using roulette selection.
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
                # Do stuff with batch. All the molecules in the batch
                # belong to the top 5 batches of size 3. The batch
                # was selected using roulette selection.
                pass
    """

    def __init__(self, filter: Selector[T], selector: Selector[T]) -> None:
        """
        Parameters:
            filter:
                Selects molecules which can be yielded by `selector`.

            selector:
                Selects batches of molecules. The batches can only
                contain molecules yielded by `filter`.
        """
        self._filter = filter
        self._selector = selector

    def select(
        self,
        population: dict[T, float],
        included_batches: "IncludedBatches" = None,
        excluded_batches: "ExcludedBatches" = None,
    ) -> Iterator[Batch[T]]:
        valid_batches = self._filter.select(
            population=population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
        valid_population = {
            record: population[record]
            for batch in valid_batches
            for record in batch
        }
        yield from self._selector.select(
            population=valid_population,
            included_batches=included_batches,
            excluded_batches=excluded_batches,
        )
