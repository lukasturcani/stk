import itertools
import typing
from collections.abc import Callable, Sequence

import numpy as np

from stk._internal.key_makers.inchi import Inchi
from stk._internal.key_makers.molecule import MoleculeKeyMaker

from .selector import Selector

T = typing.TypeVar("T")


class AboveAverage(Selector):
    """
    Yields above average batches of molecules.

    Examples:

        *Yielding Single Molecule Batches*

        Yielding molecules one at a time. For example, if molecules need
        to be selected for mutation or the next generation

        .. testcode:: yielding-single-molecule-batches

            import stk

            # Make the selector.
            above_avg = stk.AboveAverage()

            population = (
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
                ).with_fitness_value(1),
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
                ).with_fitness_value(2)
            )

            # Select the molecules.
            for selected, in above_avg.select(population):
                # Do stuff with each selected molecule.
                pass

        *Yielding Batches Holding Multiple Molecules*

        Yielding multiple molecules at once. For example, if molecules need
        to be selected for crossover.

        .. testcode:: yielding-batches-holding-multiple-molecules

            import stk

            # Make the selector.
            above_avg = stk.AboveAverage(batch_size=2)

            population = (
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
                ).with_fitness_value(1),
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
                ).with_fitness_value(2)
            )

            # Select the molecules.
            for selected1, selected2 in above_avg.select(population):
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
        fitness_modifier: Callable[[Sequence[T]], dict[T, float]]
        | None = None,
    ) -> None:
        """
        Parameters:
            num_batches:
                The number of batches to yield. If ``None`` then yielding
                will continue forever or until the generator is exhausted,
                whichever comes first.

            batch_size:
                The number of molecule records in each yielded
                :class:`.Batch`.

            duplicate_molecules:
                If ``True`` the same molecule can be yielded in more
                than one batch.

            duplicate_batches:
                If ``True`` the same batch can be yielded more than once.

            key_maker:
                Used to get the keys of molecules. If two molecules have
                the same key, they are considered duplicates.

            fitness_modifier:
                Takes the `population` on which :meth:`~Selector.select`
                is called and returns a :class:`dict`, which maps records
                in the `population` to the fitness values the
                :class:`.Selector` should use. If ``None``, the regular
                fitness values of the records are used.
        """
        if fitness_modifier is None:
            fitness_modifier = self._get_fitness_values

        super().__init__(key_maker, fitness_modifier, batch_size)

        self._duplicate_molecules = duplicate_molecules
        self._duplicate_batches = duplicate_batches
        self._num_batches = num_batches

    def _select_from_batches(self, batches, yielded_batches):
        mean = np.mean([batch.get_fitness_value() for batch in batches])
        # Yield highest fitness batches first.
        batches = sorted(batches, reverse=True)
        # Yield only batches with a fitness larger than the mean.
        batches = itertools.takewhile(
            lambda batch: batch.get_fitness_value() > mean, batches
        )
        # Yield batches which are multiple times better than the mean
        # multiple times.
        batches = (
            batch
            for batch in batches
            for _ in range(self._get_num_duplicates(batch, mean))
        )
        # If duplicate molecules are not allowed, allow only
        # batches with no yielded molecules.
        if not self._duplicate_molecules:
            batches = filter(
                yielded_batches.has_no_yielded_molecules,
                batches,
            )
        # If duplicate batches are not allowed, allow only
        # unyielded batches.
        if not self._duplicate_batches:
            batches = filter(
                yielded_batches.is_unyielded_batch,
                batches,
            )
        # Limit the number of yielded batches to _num_batches.
        yield from itertools.islice(batches, self._num_batches)

    def _get_num_duplicates(self, batch, mean):
        if self._duplicate_batches and self._duplicate_molecules:
            return int(batch.get_fitness_value() // mean)
        return 1
