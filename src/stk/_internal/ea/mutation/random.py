import typing
from collections.abc import Iterable

import numpy as np

from stk._internal.ea.mutation.mutator import MoleculeMutator
from stk._internal.ea.mutation.record import MutationRecord

T = typing.TypeVar("T")


class RandomMutator:
    """
    Use some other mutator at random.

    Examples:
        *Use One of Several Mutators at Random*

        .. testcode:: use-one-of-several-mutators-at-random

            import stk

            random_seed = 12
            mutator = stk.RandomMutator(
                mutators=(
                    stk.RandomBuildingBlock(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles='BrCCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                            stk.BuildingBlock(
                                smiles='BrCCCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                        ),
                        is_replaceable=lambda building_block: True,
                        random_seed=random_seed,
                    ),
                    stk.SimilarBuildingBlock(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles='BrCCCCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                            stk.BuildingBlock(
                                smiles='BrCCCCCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                        ),
                        is_replaceable=lambda building_block: True,
                        random_seed=random_seed,
                    ),
                    stk.RandomBuildingBlock(
                        building_blocks=(
                            stk.BuildingBlock(
                                smiles='BrCCNCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                            stk.BuildingBlock(
                                smiles='BrCNCBr',
                                functional_groups=[stk.BromoFactory()],
                            ),
                        ),
                        is_replaceable=lambda building_block: True,
                        random_seed=random_seed,
                    ),
                ),
            )
            building_block = stk.BuildingBlock(
                smiles='BrCNNCBr',
                functional_groups=[stk.BromoFactory()],
            )
            record = stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(building_block, ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            )
            # Use one of the component mutators at random.
            mutation_record1 = mutator.mutate(record)
            # A different mutator may get selected at random the second,
            # third, etc, time.
            mutation_record2 = mutator.mutate(record)

        .. testcode:: use-one-of-several-mutators-at-random
            :hide:

            _smiles = stk.Smiles()
            assert _smiles.get_key(record.get_molecule()) == 'BrCNNCCNNCBr'

            _molecule1 = (
                mutation_record1.get_molecule_record().get_molecule()
            )
            assert _smiles.get_key(_molecule1) != 'BrCNNCNNCBr'

            _molecule2 = (
                mutation_record2.get_molecule_record().get_molecule()
            )
            assert _smiles.get_key(_molecule2) != 'BrCNNCNNCBr'
    """

    def __init__(
        self,
        mutators: Iterable[MoleculeMutator[T]],
        weights: Iterable[float] | None = None,
        random_seed: int | np.random.Generator | None = None,
    ) -> None:
        """
        Parameters:
            mutators (list[MoleculeMutator[T]]):
                A selection of mutators, each time :meth:`.mutate`
                is called, one will be selected at random to perform
                the mutation.

            weights (list[float] | None):
                For each mutator, the probability that it will be chosen
                whenever :meth:`.mutate` is called.
                If ``None`` all `mutators` will have equal chance of being
                selected.

            random_seed:
                The random seed to use.
        """
        if weights is not None:
            weights = tuple(weights)

        if random_seed is None or isinstance(random_seed, int):
            random_seed = np.random.default_rng(random_seed)

        self._mutators = tuple(mutators)
        self._weights = weights
        self._generator = random_seed

    def mutate(self, record: T) -> MutationRecord[T] | None:
        """
        Return a mutant of `record`.

        Parameters:
            record:
                The molecule to be mutated.

        Returns:
            A record of the mutation or ``None``
            if `record` cannot be mutated.
        """

        mutator = self._generator.choice(
            a=self._mutators,  # type: ignore
            p=self._weights,
        )
        return mutator.mutate(record)
