"""
Random Topology Graph
=====================

"""

import numpy as np

from ....molecule_records import MoleculeRecord
from ...records import MutationRecord
from .mutator import MoleculeMutator


class RandomTopologyGraph(MoleculeMutator):
    """
    Changes topology graphs at random.

    Examples
    --------
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
        replacement_funcs,
        name='RandomTopologyGraph',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomTopology` instance.

        Parameters
        ----------
        replacement_funcs : :class:`tuple` of :class:`callable`
            Each :class:`callable` takes a single parameter, a
            :class:`.TopologyGraph`, and returns the
            :class:`.TopologyGraph` which should replace it.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """

        self._replacement_funcs = replacement_funcs
        self._name = name
        self._generator = np.random.RandomState(random_seed)

    def mutate(self, record):
        replacement_func = self._generator.choice(
            a=self._replacement_funcs,
        )
        replacement = replacement_func(record.get_topology_graph())
        return MutationRecord(
            molecule_record=MoleculeRecord(replacement),
            mutator_name=self._name,
        )
