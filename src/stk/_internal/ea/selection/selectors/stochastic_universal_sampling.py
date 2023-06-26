import typing
from collections.abc import Callable, Sequence

import numpy as np

from stk._internal.key_makers.inchi import Inchi
from stk._internal.key_makers.molecule import MoleculeKeyMaker

from .selector import Selector

T = typing.TypeVar("T")


class StochasticUniversalSampling(Selector[T]):
    """
    Yields batches of molecules through stochastic universal sampling.

    Stochastic universal sampling lays out batches along a line, with
    each batch taking up length proportional to its fitness. It
    then creates a set of evenly spaced pointers to different points
    on the line, each of which is occupied by a batch. Batches which
    are pointed to are yielded.

    This approach means weaker members of the population
    are given a greater chance to be chosen than in
    :class:`.Roulette` selection [#]_.

    References:

    .. [#] https://en.wikipedia.org/wiki/Stochastic_universal_sampling

    Examples:

        *Yielding Single Molecule Batches*

        Yielding molecules one at a time. For example, if molecules need
        to be selected for mutation or the next generation.

        .. testcode:: yielding-single-molecule-batches

            import stk

            # Make the selector.
            stochastic_sampling = stk.StochasticUniversalSampling(5)

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

            # Select the molecules.
            for selected, in stochastic_sampling.select(population):
                # Do stuff with each selected molecule.
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
        random_seed: int | np.random.Generator | None = None,
    ):
        """
        Parameters:

            num_batches:
                The number of batches to yield. If ``None`` then yielding
                will continue forever or until the generator is exhausted,
                whichever comes first.

            batch_size:
                The number of molecules yielded at once.

            duplicate_molecules:
                If ``True`` the same molecule can be yielded in more than
                one batch.

            duplicate_batches:
                If ``True`` the same batch can be yielded more than once.

            key_maker:
                Used to get the keys of molecules. If two molecules have
                the same key, they are considered duplicates.

            fitness_modifier:
                Takes the `population` on which :meth:`~.Selector.select`
                is called and returns a :class:`dict`, which maps records
                in the `population` to the fitness values the
                :class:`.Selector` should use. If ``None``, the regular
                fitness values of the records are used.

            random_seed:
                The random seed to use.
        """
        if fitness_modifier is None:
            fitness_modifier = self._get_fitness_values

        super().__init__(key_maker, fitness_modifier, batch_size)

        if random_seed is None or isinstance(random_seed, int):
            random_seed = np.random.default_rng(random_seed)

        self._generator = random_seed
        self._duplicate_molecules = duplicate_molecules
        self._duplicate_batches = duplicate_batches
        self._num_batches = (
            float("inf") if num_batches is None else num_batches
        )

    def _select_from_batches(self, batches, yielded_batches):
        batches = sorted(batches, reverse=True)

        # SUS may need to run multiple rounds if duplicate_molecules or
        # duplicate_batches is False. This is because in each round
        # you can generate multiple pointers to the same batch or to
        # batches sharing molecules. If this happens the lower fitness
        # batch will not be yielded. Instead a second round of SUS will
        # occur with any ineligible batches removed and a reduced
        # number of pointers, to account for batches yielded in the
        # previous rounds. This will repeat until the desired number
        # of batches has been yielded, or there are no more valid
        # batches.
        while batches and yielded_batches.get_num() < self._num_batches:
            yield from self._select_with_stochastic_universal_sampling(
                batches=batches,
                yielded_batches=yielded_batches,
            )

            if yielded_batches.get_num() < self._num_batches:
                if not self._duplicate_molecules:
                    batches = filter(
                        yielded_batches.has_no_yielded_molecules,
                        batches,
                    )
                if not self._duplicate_batches:
                    batches = filter(
                        yielded_batches.is_unyielded_batch,
                        batches,
                    )
                if (
                    not self._duplicate_molecules
                    or not self._duplicate_batches
                ):
                    batches = tuple(batches)

    def _select_with_stochastic_universal_sampling(
        self,
        batches,
        yielded_batches,
    ):
        total = sum(batch.get_fitness_value() for batch in batches)
        batch_positions = []
        batch_position = 0
        for batch in batches:
            batch_position += batch.get_fitness_value() / total
            batch_positions.append(batch_position)

        num_batches = min(
            self._num_batches - yielded_batches.get_num(), len(batches)
        )
        pointer_distance = 1 / num_batches
        pointers = []
        pointer = self._generator.uniform(0, pointer_distance)
        for i in range(num_batches):
            pointers.append(pointer)
            pointer += pointer_distance

        batch_index = 0
        for pointer in pointers:
            while pointer > batch_positions[batch_index]:
                batch_index += 1

            batch = batches[batch_index]

            if (
                not self._duplicate_molecules
                and yielded_batches.has_yielded_molecules(batch)
            ):
                continue

            if (
                not self._duplicate_batches
                and yielded_batches.is_yielded_batch(batch)
            ):
                continue

            yield batch
