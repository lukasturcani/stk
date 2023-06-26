from __future__ import annotations

import typing
import warnings
from collections import abc

import numpy as np

from stk._internal.building_block import BuildingBlock
from stk._internal.optimizers.null import (
    NullOptimizer,
)
from stk._internal.optimizers.optimizer import (
    Optimizer,
)
from stk._internal.reaction_factories.generic_reaction_factory import (
    GenericReactionFactory,
)
from stk._internal.reaction_factories.reaction_factory import ReactionFactory
from stk._internal.topology_graphs.edge import (
    Edge,
)
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)
from stk._internal.topology_graphs.vertex import (
    Vertex,
)

from .vertices import CycleVertex


class Macrocycle(TopologyGraph):
    """
    Represents a macrocycle topology graph.

    Building blocks with two functional groups are required
    for this topology.

    Examples:

        *Construction*

        This topology graph essentially makes a polymer chain and joins
        the ends, hence the constructor parameters allows you to
        specify the chain

        .. testcode:: construction

            import stk

            macrocycle = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCNCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=5,
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            macrocycle = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCNCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=5,
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        macrocycle.get_atoms(),
                        macrocycle.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in macrocycle.get_bonds()
                ),
            )

        *Suggested Optimization*

        For :class:`.Macrocycle` topologies, it is recommended to use
        the :class:`.MCHammer` optimizer.

        .. testcode:: suggested-optimization

            import stk

            macrocycle = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCNCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=5,
                    optimizer=stk.MCHammer(),
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            macrocycle = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCNCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='AB',
                    num_repeating_units=5,
                    optimizer=stk.MCHammer(),
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        macrocycle.get_atoms(),
                        macrocycle.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in macrocycle.get_bonds()
                ),
            )

        *Defining the Orientation of Each Building Block*

        The `orientations` parameter allows the direction of each
        building block along to the chain to be flipped

        .. testcode:: defining-the-orientation-of-each-building-block

            import stk

            bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
            bb2 = stk.BuildingBlock('BrCOCCBr', [stk.BromoFactory()])

            c1 = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(bb1, bb2),
                    repeating_unit='AB',
                    num_repeating_units=5,
                    orientations=(1, 0.5),
                ),
            )

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
            bb2 = stk.BuildingBlock('BrCOCCBr', [stk.BromoFactory()])

            macrocycle = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(bb1, bb2),
                    repeating_unit='AB',
                    num_repeating_units=5,
                    orientations=(1, 0.5),
                    random_seed=1,
                ),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        macrocycle.get_atoms(),
                        macrocycle.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in macrocycle.get_bonds()
                ),
            )

        In the above example, ``bb1`` is guaranteed to be flipped,
        ``bb2`` has a 50% chance of being flipped, each time it is
        placed on a node.

        Note that whether a building block will be flipped or not
        is decided during the initialization of :class:`.Macrocycle`

        .. testcode:: defining-the-orientation-of-each-building-block

            # cycle will always construct the same macrocycle.
            cycle = stk.macrocycle.Macrocycle(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(0.65, 0.45),
            )
            # c2 and c3 are guaranteed to be the same as they used the
            # same topology graph.
            c2 = stk.ConstructedMolecule(cycle)
            c3 = stk.ConstructedMolecule(cycle)

            # cycle2 may lead to a different polymer than chain,
            # despite being initialized with the same parameters.
            cycle2 = stk.macrocycle.Macrocycle(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(0.65, 0.45)
            )

            # c4 and c5 are guaranteed to be the same because they used
            # the same topology graph. However, they may be different
            # to c2 and c3.
            c4 = stk.ConstructedMolecule(cycle2)
            c5 = stk.ConstructedMolecule(cycle2)

        The `random_seed` parameter can be used to get reproducible
        results

        .. testcode:: defining-the-orientation-of-each-building-block

            # c6 and c7 are guaranteed to be the same, because cycle3
            # and cycle4 used the same random seed.

            cycle3 = stk.macrocycle.Macrocycle(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(0.65, 0.45),
                random_seed=4,
            )
            c6 = stk.ConstructedMolecule(cycle3)

            cycle4 = stk.macrocycle.Macrocycle(
                building_blocks=(bb1, bb2),
                repeating_unit='AB',
                num_repeating_units=5,
                orientations=(0.65, 0.45),
                random_seed=4,
            )
            c7 = stk.ConstructedMolecule(cycle4)

        *Using Numbers to Define the Repeating Unit*

        The repeating unit can also be specified through the indices of
        the building blocks

        .. testcode:: using-numbers-to-define-the-repeating-unit

            import stk

            bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
            bb2 = stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()])
            bb3 = stk.BuildingBlock('BrCNNCBr', [stk.BromoFactory()])

            # c1 and c2 are different ways to write the same thing.
            c1 = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(bb1, bb2, bb3),
                    repeating_unit='ACB',
                    num_repeating_units=3
                )
            )
            c2 = stk.ConstructedMolecule(
                topology_graph=stk.macrocycle.Macrocycle(
                    building_blocks=(bb1, bb2, bb3),
                    repeating_unit=(0, 2, 1),
                    num_repeating_units=3,
                )
            )

    """

    def __init__(
        self,
        building_blocks: tuple[BuildingBlock, ...],
        repeating_unit: typing.Union[str, tuple[int, ...]],
        num_repeating_units: int,
        orientations: tuple[float, ...] | None = None,
        random_seed: int | np.random.Generator | None = None,
        reaction_factory: ReactionFactory = GenericReactionFactory(),
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
    ) -> None:
        """
        Initialize a :class:`Macrocycle` instance.

        Parameters:
            building_blocks:
                The building blocks of the macrocycle.

            repeating_unit:
                A string specifying the repeating unit of the
                macrocycle. For example, ``'AB'`` or ``'ABB'``. The
                first building block passed to `building_blocks` is
                ``'A'`` and so on.

                The repeating unit can also be specified by the
                indices of `building_blocks`, for example ``'ABB'``
                can be written as ``(0, 1, 1)``.

            num_repeating_units:
                The number of repeating units which are used to make
                the macrocycle.

            orientations:
                For each character in the repeating unit, a value
                between ``0`` and ``1`` (both inclusive) must be given
                in a :class:`tuple`. It indicates the probability that
                each monomer will have its orientation along the chain
                flipped. If ``0`` then the monomer is guaranteed not
                to flip. If ``1`` it is guaranteed to flip. This
                allows the user to create head-to-head or head-to-tail
                chains, as well as chain with a preference for
                head-to-head or head-to-tail if a number between
                ``0`` and ``1`` is chosen. If ``None``, then ``0`` is
                picked in every case.

                It is also possible to supply an orientation for every
                vertex in the final topology graph. In this case, the
                length of `orientations` must be equal to
                ``len(repeating_unit)*num_repeating_units``.

            random_seed:
                The random seed to use when choosing random
                orientations.

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

        if orientations is None:
            orientations = tuple(0.0 for i in range(len(repeating_unit)))

        if len(orientations) == len(repeating_unit):
            orientations = orientations * num_repeating_units

        chain_length = len(repeating_unit) * num_repeating_units
        if chain_length == 2:
            warnings.warn(
                "The orientation of macrocycles with chain length "
                f"{chain_length} is not expected to provide robust "
                "alignment and bonding."
            )

        if len(orientations) != chain_length:
            raise ValueError(
                "The length of orientations must match either "
                "the length of repeating_unit or the "
                "total number of vertices."
            )

        if random_seed is None or isinstance(random_seed, int):
            random_seed = np.random.default_rng(random_seed)

        # Keep these for __repr__.
        self._repeating_unit = self._normalize_repeating_unit(
            repeating_unit=repeating_unit
        )
        self._num_repeating_units = num_repeating_units

        # Each monomer in the macrocycle is separated by angle_diff.
        angle_diff = (2 * np.pi) / chain_length
        vertices: list[CycleVertex] = []
        edges: list[Edge] = []
        choices = [True, False]
        for vertex_id, flip_chance in enumerate(orientations):
            theta = vertex_id * angle_diff
            vertices.append(
                CycleVertex(
                    id=vertex_id,
                    position=(np.cos(theta), np.sin(theta), 0),
                    flip=random_seed.choice(
                        choices,
                        p=[flip_chance, 1 - flip_chance],
                    ),
                    angle=theta,
                )
            )

            if vertex_id > 0:
                edges.append(
                    Edge(
                        id=len(edges),
                        vertex1=vertices[vertex_id - 1],
                        vertex2=vertices[vertex_id],
                    )
                )

        # Save the chosen orientations for __repr__.
        self._orientations = tuple(
            int(vertex.get_flip()) for vertex in vertices
        )

        edges.append(Edge(len(edges), vertices[0], vertices[-1]))
        super().__init__(
            building_block_vertices=self._get_building_block_vertices(
                building_blocks=building_blocks,
                vertices=vertices,
            ),
            edges=tuple(edges),
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            optimizer=optimizer,
            edge_groups=None,
        )

    def clone(self) -> Macrocycle:
        clone = self._clone()
        clone._repeating_unit = self._repeating_unit
        clone._num_repeating_units = self._num_repeating_units
        clone._orientations = self._orientations
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
        building_blocks: tuple[BuildingBlock, ...],
        vertices: abc.Sequence[CycleVertex],
    ) -> dict[BuildingBlock, abc.Sequence[Vertex]]:
        polymer = self._repeating_unit * self._num_repeating_units
        building_block_vertices: dict[BuildingBlock, list[Vertex]]
        building_block_vertices = {}
        for bb_index, vertex in zip(polymer, vertices):
            bb = building_blocks[bb_index]
            building_block_vertices[bb] = building_block_vertices.get(bb, [])
            building_block_vertices[bb].append(vertex)
        return typing.cast(
            dict[BuildingBlock, abc.Sequence[Vertex]],
            building_block_vertices,
        )

    def _get_scale(
        self,
        building_block_vertices: dict[BuildingBlock, abc.Sequence[Vertex]],
    ) -> float:
        length = len(self._repeating_unit) * self._num_repeating_units
        return (
            length
            * 0.25
            * max(bb.get_maximum_diameter() for bb in building_block_vertices)
        )

    def with_building_blocks(
        self,
        building_block_map: dict[BuildingBlock, BuildingBlock],
    ) -> Macrocycle:
        return super()._with_building_blocks(building_block_map)

    def __repr__(self) -> str:
        return (
            f"macrocycle.Macrocycle({self._repeating_unit!r}, "
            f"{self._num_repeating_units!r}, "
            f"{self._orientations!r})"
        )
