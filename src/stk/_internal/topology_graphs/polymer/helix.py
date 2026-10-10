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
from stk._internal.reaction_factories.reaction_factory import ReactionFactory
from stk._internal.topology_graphs.edge import Edge
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)
from stk._internal.topology_graphs.vertex import Vertex

from .vertices import HelixVertex, UnaligningVertex


class Helix(TopologyGraph):
    """
    Represents a linear polymer topology graph with a helical geometry.

    The model is developed from
    `Polles et. al. <https://www.nature.com/articles/ncomms7423>`_

    Examples:

        *Construction*

        Much of the same functionality as the :class:`stk.polymer.Linear`
        polymer graph is present. Helixes require building blocks with two
        functional groups.

        .. testcode:: construction

            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2],
                    repeating_unit='AB',
                    num_repeating_units=4,
                    half_rotations=2,
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2],
                    repeating_unit='AB',
                    num_repeating_units=4,
                    half_rotations=2,
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        polymer.get_atoms(),
                        polymer.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in polymer.get_bonds()
                ),
            )


        Increasing ``half_rotations`` results in a tigher helix, also note the
        change in direction with ``mirror``.

        .. testcode:: construction

            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2],
                    repeating_unit='AB',
                    num_repeating_units=4,
                    half_rotations=4,
                    mirror=True,
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2],
                    repeating_unit='AB',
                    num_repeating_units=4,
                    half_rotations=4,
                    mirror=True,
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        polymer.get_atoms(),
                        polymer.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in polymer.get_bonds()
                ),
            )

        *Suggested Optimization*

        For :class:`.Helix` topologies, it is recommend to use the
        :class:`.MCHammer` optimizer.

        .. testcode:: suggested-optimization

            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2],
                    repeating_unit='AB',
                    num_repeating_units=4,
                    half_rotations=2,
                    optimizer=stk.MCHammer(),
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2],
                    repeating_unit='AB',
                    num_repeating_units=4,
                    half_rotations=2,
                    optimizer=stk.MCHammer(),
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        polymer.get_atoms(),
                        polymer.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in polymer.get_bonds()
                ),
            )


        *Construction with Capping Units*

        Building blocks with a single functional group can
        also be provided as capping units

        .. testcode:: construction-with-capping-units

            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            bb3 = stk.BuildingBlock('BrCCN', stk.PrimaryAminoFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2, bb3],
                    repeating_unit='ABABC',
                    num_repeating_units=1,
                    half_rotations=2,
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            bb1 = stk.BuildingBlock('Nc1cccc(N)c1', stk.PrimaryAminoFactory())
            bb2 = stk.BuildingBlock('O=CC1=CC=CC=C1C=O', stk.AldehydeFactory())
            bb3 = stk.BuildingBlock('BrCCN', stk.PrimaryAminoFactory())
            polymer = stk.ConstructedMolecule(
                topology_graph=stk.polymer.Helix(
                    building_blocks=[bb1, bb2, bb3],
                    repeating_unit='ABABC',
                    num_repeating_units=1,
                    half_rotations=2,
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        polymer.get_atoms(),
                        polymer.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in polymer.get_bonds()
                ),
            )

    """

    def __init__(
        self,
        building_blocks: abc.Iterable[BuildingBlock],
        repeating_unit: str | abc.Iterable[int],
        num_repeating_units: int,
        half_rotations: int,
        mirror: bool = False,
        orientations: abc.Iterable[float] | None = None,
        random_seed: int | np.random.Generator | None = None,
        reaction_factory: ReactionFactory = GenericReactionFactory(),  # noqa: B008
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),  # noqa: B008
        scale_multiplier: float = 1.0,
    ) -> None:
        """
        Parameters:

            building_blocks (list[BuildingBlock]):
                The building blocks of the polymer.

            repeating_unit (str | list[int]):
                A string specifying the repeating unit of the polymer.
                For example, ``'AB'`` or ``'ABB'``. The first building
                block passed to `building_blocks` is ``'A'`` and so on.

                The repeating unit can also be specified by the
                indices of `building_blocks`, for example ``'ABB'``
                can be written as ``[0, 1, 1]``.

            num_repeating_units:
                The number of repeating units which are used to make
                the polymer.

            half_rotations:
                Number of rotations in the helix. The height of the helix is
                not impacted by this parameter, but comes from the size of the
                building blocks.

            mirror:
                ``True`` to form a - (counter-clockwise) twist, defaults to
                ``False`` or a + twist.

            orientations (list[float]):
                For each character in the repeating unit, a value
                between ``0`` and ``1`` (both inclusive). It indicates
                the probability that each monomer will have its
                orientation along the chain flipped. If ``0`` then the
                monomer is guaranteed not to flip. If ``1`` it is
                guaranteed to flip. This allows the user to create
                head-to-head or head-to-tail chains, as well as chain
                with a preference for head-to-head or head-to-tail if a
                number between ``0`` and ``1`` is chosen. If ``None``
                then ``0`` is picked in all cases.

                It is also possible to supply an orientation for every
                vertex in the final topology graph. In this case, the
                length of `orientations` must be equal to
                ``len(repeating_unit)*num_repeating_units``.

                If there is only one building block in the constructed
                polymer i.e. the `repeating_unit` has a length of 1 and
                `num_repeating_units` is 1, the building block will not
                be re-oriented, even if you provide a value to
                `orientations`.

            random_seed:
                The random seed to use when choosing random
                orientations.

            reaction_factory:
                The factory to use for creating reactions between
                functional groups of building blocks.

            num_processes:
                The number of parallel processes to create during
                :meth:`construct`.

            optimizer:
                Used to optimize the structure of the constructed
                molecule.

            scale_multiplier:
                Scales the positions of the vertices.

        Raises:

            :class:`ValueError`
                If the length of `orientations` is not equal in length
                to `repeating_unit` or to the total number of vertices.

        """

        self._repr = (
            f"Helix({building_blocks!r}, "
            f"{repeating_unit!r}, {num_repeating_units!r}, {half_rotations!r},"
            f" {mirror!r})"
        )

        if not isinstance(repeating_unit, str):
            repeating_unit = tuple(repeating_unit)

        if not isinstance(repeating_unit, str):
            repeating_unit = tuple(repeating_unit)

        if orientations is None:
            orientations = tuple(0.0 for _ in range(len(repeating_unit)))
        else:
            orientations = tuple(orientations)

        if len(orientations) == len(repeating_unit):
            orientations = orientations * num_repeating_units

        polymer_length = len(repeating_unit) * num_repeating_units
        if len(orientations) != polymer_length:
            raise ValueError(
                "The length of orientations must match either "
                "the length of repeating_unit or the "
                "total number of vertices."
            )

        # Keep these for __repr__.
        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit
        )
        self._num_repeating_units = num_repeating_units

        try:
            head, *body, tail = orientations
            vertices_and_edges = self._get_vertices_and_edges(
                head_orientation=head,
                body_orientations=body,
                tail_orientation=tail,
                random_seed=random_seed,
                mirror=-1 if mirror else 1,
                half_rotations=half_rotations,
            )
            vertices = vertices_and_edges.vertices
            edges = vertices_and_edges.edges

        except ValueError:
            vertices = (UnaligningVertex(0, (0.0, 0.0, 0.0), False),)
            edges = ()

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(int(v.get_flip()) for v in vertices)

        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                building_blocks=tuple(building_blocks),
                vertices=vertices,
            ),
            edges=edges,
            reaction_factory=reaction_factory,
            construction_stages=(),
            optimizer=optimizer,
            num_processes=num_processes,
            scale_multiplier=scale_multiplier,
        )

    @staticmethod
    def _get_vertices_and_edges(
        head_orientation: float,
        body_orientations: abc.Iterable[float],
        tail_orientation: float,
        random_seed: int | np.random.Generator | None,
        mirror: int,
        half_rotations: int,
    ) -> "_VerticesAndEdges":
        """
        Get the vertices and edges of the topology graph.

        Parameters:

            head_orientation:
                The probability that the head vertex will do flipping.

            body_orientations:
                For each body vertex, the probability that it will do
                flipping.

            tail_orientation:
                The probability that the tail vertex will do flipping.

            random_seed:
                The random seed to use.

            mirror:
                -1 if the mirror structure is desired.

            half_rotations:
                Number of rotations in helix. The heigh of the helix is
                unchanged.

        Returns:

            The vertices and edges of the topology graph.

        """

        max_angle = np.pi * half_rotations
        radius = 1
        height = radius * 1
        num_bbs = len(list(body_orientations)) + 1
        v_z = height / max_angle  # z(t) = v_z * t
        bb_distance = 1

        # approximate angle step
        half_step_hsp = np.sqrt(bb_distance**2 / (v_z**2 + radius**2))
        step_hsp = (max_angle - 2 * half_step_hsp) / (num_bbs - 1)

        if random_seed is None or isinstance(random_seed, int):
            random_seed = np.random.default_rng(random_seed)

        choices = [True, False]
        vertices: list[HelixVertex | UnaligningVertex] = [
            HelixVertex(
                id=0,
                position=np.array([radius, 0, 0 * mirror]),
                flip=bool(
                    random_seed.choice(
                        a=choices,
                        p=[head_orientation, 1 - head_orientation],
                    )
                ),
            ),
        ]
        edges: list[Edge] = []
        for i, p in enumerate(body_orientations, 1):
            curr_angle = half_step_hsp + i * step_hsp
            flip = bool(random_seed.choice(choices, p=[p, 1 - p]))
            vertices.append(
                HelixVertex(
                    i,
                    np.array(
                        [
                            radius * np.cos(curr_angle),
                            radius * np.sin(curr_angle),
                            curr_angle * v_z * mirror,
                        ]
                    ),
                    flip,
                ),
            )
            edges.append(Edge(len(edges), vertices[i - 1], vertices[i]))

        vertices.append(
            HelixVertex(
                id=len(vertices),
                position=np.array(
                    [
                        radius * np.cos(max_angle),
                        radius * np.sin(max_angle),
                        height * mirror,
                    ]
                ),
                flip=bool(
                    random_seed.choice(
                        a=choices,
                        p=[tail_orientation, 1 - tail_orientation],
                    )
                ),
            ),
        )

        edges.append(Edge(len(edges), vertices[-2], vertices[-1]))

        return _VerticesAndEdges(
            vertices=tuple(vertices),
            edges=tuple(edges),
        )

    def clone(self) -> typing.Self:
        clone = self._clone()
        clone._repr = self._repr
        clone._repeating_unit = self._repeating_unit
        clone._num_repeating_units = self._num_repeating_units
        clone._orientations = self._orientations
        return clone

    @staticmethod
    def _normalize_repeating_unit(
        repeating_unit: str | tuple[int, ...],
    ) -> tuple[int, ...]:
        if isinstance(repeating_unit, tuple):
            return repeating_unit

        base = ord("A")
        return tuple(ord(letter) - base for letter in repeating_unit)

    def _get_building_block_vertices(
        self,
        building_blocks: tuple[BuildingBlock, ...],
        vertices: tuple[HelixVertex | UnaligningVertex, ...],
    ) -> dict[BuildingBlock, abc.Sequence[Vertex]]:
        polymer = self._repeating_unit * self._num_repeating_units

        building_block_vertices: dict[
            BuildingBlock, list[HelixVertex | UnaligningVertex]
        ] = {}

        for bb_index, vertex in zip(polymer, vertices):
            bb = building_blocks[bb_index]
            building_block_vertices[bb] = building_block_vertices.get(bb, [])
            building_block_vertices[bb].append(vertex)

        return self._with_unaligning_vertices(
            building_block_vertices=building_block_vertices,
        )

    @staticmethod
    def _with_unaligning_vertices(
        building_block_vertices: dict[
            BuildingBlock, list[HelixVertex | UnaligningVertex]
        ],
    ) -> dict[BuildingBlock, abc.Sequence[Vertex]]:
        clone: dict[BuildingBlock, abc.Sequence[Vertex]]
        clone = {}
        terminal_ids = {
            0,
            max(
                vertex.get_id()
                for vertex_list in building_block_vertices.values()
                for vertex in vertex_list
            ),
        }
        for (
            building_block,
            vertices,
        ) in building_block_vertices.items():
            # Building blocks with 1 placer, cannot be aligned on
            # linear vertices and must therefore use an
            # UnaligningVertex. Building blocks with 1 placer can be
            # placed on terminal vertices (HeadVertex or TailVertex).
            # This can be discerned based on the knowledge that the
            # first and last vertex are the Head and Tail,
            # respectively.
            if building_block.get_num_placers() == 1:
                clone[building_block] = tuple(
                    (
                        UnaligningVertex(
                            id=vertex.get_id(),
                            position=vertex.get_position(),
                            flip=vertex.get_flip(),
                        )
                        if vertex.get_id() not in terminal_ids
                        else vertex
                    )
                    for vertex in vertices
                )
            else:
                clone[building_block] = vertices

        return clone

    @staticmethod
    def _get_scale(
        building_block_vertices: dict[BuildingBlock, abc.Sequence[Vertex]],
        scale_multiplier: float,
    ) -> float:
        return scale_multiplier * max(
            bb.get_maximum_diameter() for bb in building_block_vertices
        )

    def with_building_blocks(
        self,
        building_block_map: dict[BuildingBlock, BuildingBlock],
    ) -> typing.Self:
        return self.clone()._with_building_blocks(building_block_map)

    def __repr__(self) -> str:
        return self._repr


@dataclass(frozen=True)
class _VerticesAndEdges:
    vertices: tuple[HelixVertex | UnaligningVertex, ...]
    edges: tuple[Edge, ...]
