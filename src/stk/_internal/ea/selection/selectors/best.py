import itertools
import typing
from collections.abc import Callable, Iterable, Iterator, Sequence

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import Batch
from stk._internal.ea.selection.selectors.yielded_batches import YieldedBatches
from stk._internal.key_makers.inchi import Inchi
from stk._internal.key_makers.molecule import MoleculeKeyMaker

from .selector import Selector

T = typing.TypeVar("T", bound=MoleculeRecord)


class Best(Selector[T]):
    """
    Selects batches of molecules, highest fitness value first.

    Examples:
        *Yielding Single Molecule Batches*

        Yielding molecules one at a time. For example, if molecules need
        to be selected for mutation or the next generation.

        .. testcode:: yielding-single-molecule-batches

            import stk

            # Make the selector.
            best = stk.Best()

            population = {
                stk.MoleculeRecord(
                    topology_graph=stk.polymer.Linear(
                        building_blocks=[
                            stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                        ],
                        repeating_unit='A',
                        num_repeating_units=2,
                    ),
                ): 1,
                stk.MoleculeRecord(
                    topology_graph=stk.polymer.Linear(
                        building_blocks=[
                            stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                        ],
                        repeating_unit='A',
                        num_repeating_units=2,
                    ),
                ): 2
            }

            # Select the molecules.
            for selected, in best.select(population):
                # Do stuff with each selected molecule.
                pass

        *Yielding Batches Holding Multiple Molecules*

        Yielding multiple molecules at once. For example, if molecules need
        to be selected for crossover.

        .. testcode:: yielding-batches-holding-multiple-molecules

            import stk

            # Make the selector.
            best = stk.Best(batch_size=2)

            population = {
                stk.MoleculeRecord(
                    topology_graph=stk.polymer.Linear(
                        building_blocks=[
                            stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                        ],
                        repeating_unit='A',
                        num_repeating_units=2,
                    ),
                ): 1,
                stk.MoleculeRecord(
                    topology_graph=stk.polymer.Linear(
                        building_blocks=[
                            stk.BuildingBlock('BrCCBr', stk.BromoFactory()),
                        ],
                        repeating_unit='A',
                        num_repeating_units=2,
                    ),
                ): 2
            }

            # Select the molecules.
            for selected1, selected2 in best.select(population):
                # Do stuff with the selected molecules.
                pass
    """

    def __init__(
        self,
        num_batches: int | None = None,
        batch_size: int = 1,
        duplicate_molecules: bool = True,
        duplicate_batches: bool = True,
        key_maker: MoleculeKeyMaker = Inchi(),
        fitness_modifier: Callable[
            [dict[T, float]], dict[T, float]
        ] = lambda x: x,
    ) -> None:
        """
        Parameters:
            num_batches:
                The number of batches to yield. If ``None`` then yielding
                will continue forever or until the generator is exhausted,
                whichever comes first.

            batch_size:
                The number of molecule records yielded at once.

            duplicate_molecules:
                If ``True`` the same molecule can be yielded in more than
                one batch.

            duplicate_batches:
                If ``True`` the same batch can be yielded more than once.
                Duplicate batches can occur if the same molecule is found
                multiple times in a population.

            key_maker:
                Used to get the keys of molecules. If two molecules have
                the same key, they are considered duplicates.

            fitness_modifier:
                Takes the `population` on which :meth:`~.Selector.select`
                is called and returns a :class:`dict`, which maps records
                in the `population` to the fitness values the
                :class:`.Selector` should use.
        """
        super().__init__(key_maker, fitness_modifier, batch_size)
        self._duplicate_molecules = duplicate_molecules
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches
        self._batch_size = batch_size

    def _select_from_batches(
        self,
        batches: Sequence[Batch[T]],
        yielded_batches: YieldedBatches[T],
    ) -> Iterator[Batch[T]]:
        selected_batches: Iterable[Batch[T]] = sorted(batches, reverse=True)

        if not self._duplicate_molecules:
            selected_batches = filter(
                yielded_batches.has_no_yielded_molecules,
                selected_batches,
            )

        if not self._duplicate_batches:
            selected_batches = filter(
                yielded_batches.is_unyielded_batch,
                selected_batches,
            )

        yield from itertools.islice(selected_batches, self._num_batches)
