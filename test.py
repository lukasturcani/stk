import typing
from collections import abc

import numpy as np

import stk
from stk._internal.topology_graphs.utilities import (
    _EdgeSorter,
    _FunctionalGroupSorter,
)


class SingleVertex(stk.Vertex):
    """ """

    def place_building_block(
        self,
        building_block: stk.BuildingBlock,
        edges: tuple[stk.Edge, ...],
    ) -> np.ndarray:
        return building_block.with_centroid(
            np.array([0.0, 0.0, 0.0])
        ).get_position_matrix()

    def map_functional_groups_to_edges(
        self,
        building_block: stk.BuildingBlock,
        edges: tuple[stk.Edge, ...],
    ) -> dict[int, int]:
        # The idea is to order the functional groups in building_block
        # by their angle with the vector running from the placer
        # centroid to fg0, going in the clockwise direction.
        # The edges are also ordered by their angle with the vector
        # running from the edge centroid to the aligner_edge,
        # going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg_sorter = _FunctionalGroupSorter(building_block)
        edge_sorter = _EdgeSorter(
            edges=edges,
            aligner_edge=edges[0],
            axis=fg_sorter.get_axis(),
        )
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(
                fg_sorter.get_items(),
                edge_sorter.get_items(),
            )
        }


class InternalReaction(stk.TopologyGraph):
    """ """

    def __init__(
        self,
        building_block: stk.BuildingBlock,
        num_reactions: int,
        num_processes: int = 1,
        optimizer: stk.Optimizer = stk.NullOptimizer(),
        scale_multiplier: float = 1.0,
        reaction_factory: stk.ReactionFactory = stk.GenericReactionFactory(),
    ) -> None:
        """
        Parameters:


        Raises:

            :class:`ValueError`
                If the number of reactions does not match the number of pairs
                of functional groups.

        """

        self._repr = f"Internal({building_block!r}, {num_reactions!r})"
        self._num_reactions = num_reactions
        num_functional_groups = building_block.get_num_functional_groups()
        if num_functional_groups != 2 * self._num_reactions:
            raise ValueError(
                f"The number of reactions {num_reactions} does match the number"
                f" of pairs of functional groups ({self._num_arms})."
            )

        vertices = (SingleVertex(0, (0, 0, 0)),)
        edges = []
        for potential_edge in range(self._num_reactions):
            edges.append(
                stk.Edge(
                    id=potential_edge,
                    vertex1=vertices[0],
                    vertex2=vertices[0],
                )
            )
        edges = tuple(edges)

        super().__init__(
            building_block_vertices={building_block: vertices},
            edges=edges,
            reaction_factory=reaction_factory,
            construction_stages=(),
            optimizer=optimizer,
            num_processes=num_processes,
            scale_multiplier=scale_multiplier,
        )

    def clone(self) -> typing.Self:
        clone = self._clone()
        clone._repr = self._repr
        clone._num_reactions = self._num_reactions
        return clone

    @staticmethod
    def _get_scale(
        building_block_vertices: dict[
            stk.BuildingBlock, abc.Sequence[stk.Vertex]
        ],
        scale_multiplier: float,
    ) -> float:
        return 1

    def with_building_blocks(
        self,
        building_block_map: dict[stk.BuildingBlock, stk.BuildingBlock],
    ) -> typing.Self:
        return self.clone()._with_building_blocks(building_block_map)

    def __repr__(self) -> str:
        return self._repr


cycle = stk.ConstructedMolecule(
    topology_graph=stk.macrocycle.Macrocycle(
        building_blocks=(
            stk.BuildingBlock(
                smiles="[Br]CC[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
            stk.BuildingBlock(
                smiles="[Br]C(CCI)C[Br]",
                functional_groups=[stk.BromoFactory()],
            ),
        ),
        repeating_unit="ABAAA",
        num_repeating_units=1,
    ),
)
axle = stk.ConstructedMolecule(
    topology_graph=stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()]),
            stk.BuildingBlock("BrCOCI", [stk.BromoFactory()]),
        ),
        repeating_unit="ABABC",
        num_repeating_units=1,
    )
)
rotaxane = stk.ConstructedMolecule(
    topology_graph=stk.rotaxane.NRotaxane(
        axle=stk.BuildingBlock.init_from_molecule(axle),
        cycles=(stk.BuildingBlock.init_from_molecule(cycle),),
        repeating_unit="A",
        num_repeating_units=1,
    ),
)
rotaxane.write("rotaxane.mol")


onerotax = stk.ConstructedMolecule(
    topology_graph=InternalReaction(
        building_block=stk.BuildingBlock.init_from_molecule(
            molecule=rotaxane,
            functional_groups=(stk.IodoFactory(),),
        ),
        num_reactions=1,
    ),
)
onerotax.write("1_rotaxane.mol")
