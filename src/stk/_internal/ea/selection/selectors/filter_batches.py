import typing
from collections.abc import Set

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import BatchKey

from .selector import Selector

T = typing.TypeVar("T", bound=MoleculeRecord)


class FilterBatches(Selector[T]):
    """
    Allows a :class:`.Selector` to select only some batches.

    Examples:

        *Using a Selection Algorithm on a Subset of Batches*

        You only want the :class:`.Best` 10 batches to participate in
        :class:`.Roulette`

        .. testcode:: using-a-selection-algorithm-on-a-subset-of-batches

            import stk

            # Select the 10 best batches first, and then run roulette on
            # those.
            selector = stk.FilterBatches(
                filter=stk.Best(10),
                selector=stk.Roulette(7),
            )

            population = tuple(
                stk.MoleculeRecord(
                    topology_graph=stk.polymer.Linear(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles='BrCCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                        ),
                        repeating_unit='A',
                        num_repeating_units=2,
                    ),
                ).with_fitness_value(i)
                for i in range(100)
            )

            for batch in selector.select(population):
                # Do stuff with batch. It is one of the 10 best batches and
                # was selected using roulette selection.
                pass
    """

    def __init__(self, filter: Selector[T], selector: Selector[T]) -> None:
        """
        Parameters:
            filter:
                Selects batches which can be yielded by `selector`.

            selector:
                Selects batches, but only if they were also selected by
                `filter`.
        """
        self._filter = filter
        self._selector = selector

    def select(
        self,
        population: dict[T, float],
        included_batches: Set[BatchKey] | None = None,
        excluded_batches: Set[BatchKey] | None = None,
    ):
        allowed = {
            batch.get_identity_key()
            for batch in self._filter.select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            )
        }
        yield from self._selector.select(
            population=population,
            included_batches=allowed,
            excluded_batches=excluded_batches,
        )
