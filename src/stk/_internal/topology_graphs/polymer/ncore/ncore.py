import typing
from collections import abc
from dataclasses import dataclass

import numpy as np

from stk._internal.building_block import BuildingBlock
from stk._internal.optimizers.null import NullOptimizer
from stk._internal.optimizers.optimizer import Optimizer
from stk._internal.reaction_factories.generic_reaction_factory import (
    GenericReactionFactory,
)
from stk._internal.reaction_factories.reaction_factory import (
    ReactionFactory,
)
from stk._internal.topology_graphs.edge import Edge
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)
from stk._internal.topology_graphs.vertex import Vertex

from ..vertices import CoreVertex, SubstituentVertex


class NCore(TopologyGraph):
    """
    Represents a linear polymer topology graph.

    Building blocks with two functional groups are required, unless the
    building block's position is specified to only be at the capping
    positions.

    Examples:

    """

    def __init__(
        self,
        core_building_block: BuildingBlock,
        arm_building_blocks: abc.Iterable[BuildingBlock],
        repeating_unit: str | abc.Iterable[int],
        reaction_factory: ReactionFactory = GenericReactionFactory(),
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
    ) -> None:
        """
        Parameters:

            core_building_block (list[BuildingBlock]):
                The central building block.

            arm_building_blocks (list[BuildingBlock]):
                The building blocks to place on arms.

            repeating_unit (str | list[int]):
                A string specifying the repeating unit of the polymer.
                For example, ``'AB'`` or ``'ABB'``. The first building
                block passed to `building_blocks` is ``'A'`` and so on.

                The repeating unit can also be specified by the
                indices of `building_blocks`, for example ``'ABB'``
                can be written as ``[0, 1, 1]``.

            reaction_factory:
                The factory to use for creating reactions between
                functional groups of building blocks.

            num_processes:
                The number of parallel processes to create during
                :meth:`construct`.

            optimizer:
                Used to optimize the structure of the constructed
                molecule.

        Raises:

            :class:`ValueError`
                If the length of `orientations` is not equal in length
                to `repeating_unit` or to the total number of vertices.

        """

        self._repr = (
            f"NCore({core_building_block!r}, {arm_building_blocks!r}, "
            f"{repeating_unit!r})"
        )

        if not isinstance(repeating_unit, str):
            repeating_unit = tuple(repeating_unit)

        self._num_arms = core_building_block.get_num_functional_groups()
        if self._num_arms % len(repeating_unit) != 0:
            raise ValueError(
                f"The repeating unit {repeating_unit} does not fit "
                "evenly onto the core building block with "
                f"{self._num_arms} functional groups."
            )

        # Keep these for __repr__.
        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit
        )
        self._num_repeating_units = int(self._num_arms / len(repeating_unit))

        vertices_and_edges = self._get_vertices_and_edges(self._num_arms)
        vertices = vertices_and_edges.vertices
        edges = vertices_and_edges.edges

        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                core_building_block=core_building_block,
                arm_building_blocks=tuple(arm_building_blocks),
                vertices=vertices,
            ),
            edges=edges,
            reaction_factory=reaction_factory,
            construction_stages=(),
            optimizer=optimizer,
            num_processes=num_processes,
        )

    @staticmethod
    def _get_vertices_and_edges(num_arms) -> "_VerticesAndEdges":
        """
        Get the vertices and edges of the topology graph.

        Parameters:


        Returns:

            The vertices and edges of the topology graph.

        """

        if num_arms == 1:
            arm_positions = [np.array((1, 0, 0))]
        elif num_arms == 2:
            arm_positions = [np.array((1, 0, 0)), np.array((-1, 0, 0))]
        else:
            thetas = np.linspace(0, 2 * np.pi, num_arms + 1)[:-1]
            x_points = np.cos(thetas)
            y_points = np.sin(thetas)
            arm_positions = [
                np.array((x, y, 0)) for x, y in zip(x_points, y_points)
            ]

        vertices: list[Vertex] = [CoreVertex(0, (0, 0, 0))]
        edges: list[Edge] = []
        for i, pos in enumerate(arm_positions):
            vertices.append(SubstituentVertex(i + 1, pos))
            edges.append(Edge(len(edges), vertices[0], vertices[-1]))

        return _VerticesAndEdges(
            vertices=tuple(vertices),
            edges=tuple(edges),
        )

    def clone(self) -> typing.Self:
        clone = self._clone()
        clone._repr = self._repr
        clone._repeating_unit = self._repeating_unit
        clone._num_repeating_units = self._num_repeating_units
        return clone

    @staticmethod
    def _normalize_repeating_unit(
        repeating_unit: typing.Union[str, tuple[int, ...]],
    ) -> tuple[int, ...]:
        if isinstance(repeating_unit, tuple):
            return repeating_unit

        base = ord("A")
        return tuple(ord(letter) - base for letter in repeating_unit)

    def _get_building_block_vertices(
        self,
        core_building_block: BuildingBlock,
        arm_building_blocks: tuple[BuildingBlock, ...],
        vertices: tuple[Vertex, ...],
    ) -> dict[BuildingBlock, abc.Sequence[Vertex]]:
        building_block_vertices: dict[BuildingBlock, list[Vertex]] = {}
        building_block_vertices[core_building_block] = [vertices[0]]

        bb_order = self._repeating_unit * self._num_repeating_units
        for bb_index, vertex in zip(bb_order, vertices[1:]):
            bb = arm_building_blocks[bb_index]
            building_block_vertices[bb] = building_block_vertices.get(bb, [])
            building_block_vertices[bb].append(vertex)

        # Have to do this to match the typing in the return statement
        # (Sequence), suggestions for improvement would be great.
        return {
            i: tuple(building_block_vertices[i])
            for i in building_block_vertices
        }

    def _get_scale(
        self,
        building_block_vertices: dict[BuildingBlock, abc.Sequence[Vertex]],
    ) -> float:
        return max(bb.get_maximum_diameter() for bb in building_block_vertices)

    def with_building_blocks(
        self,
        building_block_map: dict[BuildingBlock, BuildingBlock],
    ) -> typing.Self:
        return self.clone()._with_building_blocks(building_block_map)

    def __repr__(self) -> str:
        return self._repr


@dataclass(frozen=True)
class _VerticesAndEdges:
    vertices: tuple[Vertex, ...]
    edges: tuple[Edge, ...]
