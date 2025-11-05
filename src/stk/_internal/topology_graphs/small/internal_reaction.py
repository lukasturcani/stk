import typing
from collections import abc

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

from .vertices import SingleVertex


class InternalReaction(TopologyGraph):
    """Represents an intramolecular reaction(s).

    Here, we designed a way to take a single building block (maybe constructed
    from another series of building blocks) and perform any number of
    intramolecular reactions.

    Examples:

        *Construction*

        The number of functional groups in the core building block
        define the topology graph


        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            building_block = stk.BuildingBlock(
                smiles="CCCC(CCCCI)CCCCCC(I)CC",
                functional_groups=stk.IodoFactory(),
            )


            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        building_block.get_atoms(),
                        building_block.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in building_block.get_bonds()
                ),
            )


        .. testcode:: construction-internal

            import stk

            building_block = stk.BuildingBlock(
                smiles="CCCC(CCCCI)CCCCCC(I)CC",
                functional_groups=stk.IodoFactory(),
            )

            reacted = stk.ConstructedMolecule(stk.small.InternalReaction(
                building_block=building_block,
                num_reactions=1
            ))

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            building_block = stk.BuildingBlock(
                smiles="CCCC(CCCCI)CCCCCC(I)CC",
                functional_groups=stk.IodoFactory(),
            )

            reacted = stk.ConstructedMolecule(stk.small.InternalReaction(
                building_block=building_block,
                num_reactions=1
            ))

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        reacted.get_atoms(),
                        reacted.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in reacted.get_bonds()
                ),
            )

        And multiple reactions can be performed, where the closest pairs should
        react.

        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            building_block = stk.BuildingBlock(
                smiles="ICCCCC(CCCI)CCCCCC(I)CCI",
                functional_groups=stk.IodoFactory(),
            )

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        building_block.get_atoms(),
                        building_block.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in building_block.get_bonds()
                ),
            )

        .. testcode:: construction

            import stk

            building_block = stk.BuildingBlock(
                smiles="ICCCCC(CCCI)CCCCCC(I)CCI",
                functional_groups=stk.IodoFactory(),
            )
            reacted = stk.ConstructedMolecule(stk.small.InternalReaction(
                building_block=building_block,
                num_reactions=2,
            ))


        .. moldoc::

            import moldoc.molecule as molecule
            import stk

            building_block = stk.BuildingBlock(
                smiles="ICCCCC(CCCI)CCCCCC(I)CCI",
                functional_groups=stk.IodoFactory(),
            )

            reacted = stk.ConstructedMolecule(stk.small.InternalReaction(
                building_block=building_block,
                num_reactions=2
            ))

            moldoc_display_molecule = molecule.Molecule(
                atoms=(
                    molecule.Atom(
                        atomic_number=atom.get_atomic_number(),
                        position=position,
                    ) for atom, position in zip(
                        reacted.get_atoms(),
                        reacted.get_position_matrix(),
                    )
                ),
                bonds=(
                    molecule.Bond(
                        atom1_id=bond.get_atom1().get_id(),
                        atom2_id=bond.get_atom2().get_id(),
                        order=bond.get_order(),
                    ) for bond in reacted.get_bonds()
                ),
            )

        *Suggested Optimization*

        For :class:`.InternalReaction` topologies, no :mod:`stk` optimizer is
        appropriate, as they all assume distinct building blocks.

    """

    def __init__(
        self,
        building_block: BuildingBlock,
        num_reactions: int,
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
        scale_multiplier: float = 1.0,
        reaction_factory: ReactionFactory = GenericReactionFactory(),
    ) -> None:
        """
        Parameters:

            building_block:
                The only building block.

            num_reactions:
                The number of internal reactions to run.

            num_processes:
                The number of parallel processes to create during
                :meth:`construct`.

            optimizer:
                Used to optimize the structure of the constructed
                molecule.

            scale_multiplier:
                Scales the positions of the vertices.

            reaction_factory:
                The factory to use for creating reactions between
                functional groups of building blocks.

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
                f" of pairs of functional groups ({num_functional_groups})."
            )

        vertices = (SingleVertex(0, (0, 0, 0)),)
        edges = []
        for potential_edge in range(self._num_reactions):
            edges.append(
                Edge(
                    id=potential_edge,
                    vertex1=vertices[0],
                    vertex2=vertices[0],
                )
            )

        super().__init__(
            building_block_vertices={building_block: vertices},
            edges=tuple(edges),
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
        building_block_vertices: dict[BuildingBlock, abc.Sequence[Vertex]],
        scale_multiplier: float,
    ) -> float:
        return 1

    def with_building_blocks(
        self,
        building_block_map: dict[BuildingBlock, BuildingBlock],
    ) -> typing.Self:
        return self.clone()._with_building_blocks(building_block_map)

    def __repr__(self) -> str:
        return self._repr
