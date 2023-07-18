from collections.abc import Callable, Iterable
from typing import Any

import numpy as np

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.mutation.record import MutationRecord
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)


class RandomTopologyGraph:
    """
    Changes topology graphs at random.

    Examples:
        *Constructed Molecule Mutation*

        .. testcode:: constructed-molecule-mutation

            import stk

            # Create a molecule which is to be mutated.
            bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
            bb2 = stk.BuildingBlock(
                smiles='O=CCC(C=O)CC=O',
                functional_groups=[stk.AldehydeFactory()],
            )
            cage = stk.MoleculeRecord(
                topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
            )

            # Create functions which replace the topology graph.
            replacement_funcs = (
                lambda graph:
                    stk.cage.TwoPlusThree(graph.get_building_blocks()),
                lambda graph:
                    stk.cage.EightPlusTwelve(graph.get_building_blocks()),
                lambda graph:
                    stk.cage.TwentyPlusThirty(graph.get_building_blocks()),
            )

            # Create the mutator.
            random_topology = stk.RandomTopologyGraph(replacement_funcs)

            # Mutate a molecule.
            mutation_record1 = random_topology.mutate(cage)

            # Mutate the molecule a second time.
            mutation_record2 = random_topology.mutate(cage)

    """

    def __init__(
        self,
        replacement_funcs: Iterable[Callable[[TopologyGraph], TopologyGraph]],
        name: str = "RandomTopologyGraph",
        random_seed: int | np.random.Generator | None = None,
    ) -> None:
        """
        Parameters:
            replacement_funcs \
(list[collections.abc.Callable[[TopologyGraph], TopologyGraph]]):
                Functions which take a topology graph and return
                its replacement. One is selected at random each
                time :meth:`.mutate` is called.

            name:
                A name to help identify the mutator instance.

            random_seed:
                The random seed to use.

        """

        if random_seed is None or isinstance(random_seed, int):
            random_seed = np.random.default_rng(random_seed)

        self._replacement_funcs = tuple(replacement_funcs)
        self._name = name
        self._generator = random_seed

    def mutate(
        self,
        record: MoleculeRecord[Any],
    ) -> MutationRecord[MoleculeRecord[Any]]:
        """
        Return a mutant of `record`.

        Parameters:
            record:
                The molecule to be mutated.

        Returns:
            A record of the mutation or ``None``
            if `record` cannot be mutated.
        """
        replacement_func = self._generator.choice(
            a=self._replacement_funcs,  # type: ignore
        )
        replacement = replacement_func(record.get_topology_graph())
        return MutationRecord(
            molecule_record=MoleculeRecord(replacement),
            mutator_name=self._name,
        )
