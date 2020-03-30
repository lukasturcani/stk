"""
Random Topology Graph
=====================

"""

import numpy as np

from .mutator import ConstructedMoleculeMutator
from ..record import ConstructedMoleculeMutationRecord
from ....molecule_records import ConstructedMoleculeRecord


class RandomTopologyGraph(ConstructedMoleculeMutator):
    """
    Changes topology graphs at random.

    Examples
    --------
    *Constructed Molecule Mutation*

    .. code-block:: python

        import stk

        # Create a molecule which is to be mutated.
        bb1 = stk.BuildingBlock('NCCN', [stk.PrimaryAminoFactory()])
        bb2 = stk.BuildingBlock(
            smiles='O=CCC(=O)CC=O',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage = stk.ConstructedMolecule(stk.cage.FourPlusSix((bb1, bb2))

        # Create topologies used for substitution.
        topology_graphs = (
            stk.cage.TwoPlusThree(),
            stk.cage.EightPlusTwelve(),
            stk.cage.TwentyPlusThirty()
        )

        # Create the mutator.
        random_topology = stk.RandomTopologyGraph(topology_graphs)

        # Mutate a molecule.
        mutation_record1 = random_topology.mutate(cage)

        # Mutate the molecule a second time.
        mutation_record2 = random_topology.mutate(cage)

    """

    def __init__(
        self,
        topology_graphs,
        name='RandomTopologyGraph',
        random_seed=None,
    ):
        """
        Initialize a :class:`.RandomTopology` instance.

        Parameters
        ----------
        topology_graphs : :class:`tuple` of :class:`.TopologyGraph`
            Holds the topology instances from which one is
            selected at random to form a new molecule.

        name : :class:`str`, optional
            A name to help identify the mutator instance.

        random_seed : :class:`bool`, optional
            The random seed to use.

        """

        self._topology_graphs = topology_graphs
        self._name = name
        self._generator = np.random.RandomState(random_seed)

    def _mutate(self, record):
        topology_graph = self._generator.choice(self._topology_graphs)
        replacement = topology_graph.with_building_blocks({
            building_block1: building_block2
            for building_block1, building_block2
            in zip(
                topology_graph.get_building_blocks(),
                record.get_topology_graph().get_building_blocks(),
            )
        })
        return ConstructedMoleculeMutationRecord(
            molecule_record=ConstructedMoleculeRecord(replacement),
            mutator_name=self._name,
        )
