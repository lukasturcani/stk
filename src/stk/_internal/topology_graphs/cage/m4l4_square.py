from collections.abc import Iterable

import numpy as np

from stk._internal.building_block import BuildingBlock
from stk._internal.optimizers.null import NullOptimizer
from stk._internal.optimizers.optimizer import Optimizer
from stk._internal.reaction_factories.generic_reaction_factory import (
    GenericReactionFactory,
)
from stk._internal.reaction_factories.reaction_factory import ReactionFactory
from stk._internal.topology_graphs.edge import Edge

from .cage import Cage
from .vertices import AngledVertex, LinearVertex


class M4L4Square(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(2)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles=(
                'C1=CC(=CC=C1C2=CC=NC=C2)C3=CC=NC=C3'
            ),
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L4Square(
                corners=bb1,
                linkers=bb2,
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        1
                        if bond.get_order() == 9
                        else bond.get_order()
                    ),
                ) for bond in cage.get_bonds()
            ),
        )

    :class:`.MCHammer` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(2)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles=(
                'C1=CC(=CC=C1C2=CC=NC=C2)C3=CC=NC=C3'
            ),
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M4L4Square(
                corners=bb1,
                linkers=bb2,
                optimizer=stk.MCHammer(),
            ),
        )

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                ) for atom, position in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                )
            ),
            bonds=(
                molecule.Bond(
                    atom1_id=bond.get_atom1().get_id(),
                    atom2_id=bond.get_atom2().get_id(),
                    order=(
                        1
                        if bond.get_order() == 9
                        else bond.get_order()
                    ),
                ) for bond in cage.get_bonds()
            ),
        )

    Both `corner` and `linker` vertices require building blocks with
    two functional groups for this topology. This class replaces the
    `building_blocks` parameter with the `corner` and `linker`
    parameters.

    See :class:`.Cage` for more details and examples.

    """

    def __init__(
        self,
        corners: (
            Iterable[BuildingBlock] | dict[BuildingBlock, tuple[int, ...]]
        ),
        linkers: (
            Iterable[BuildingBlock] | dict[BuildingBlock, tuple[int, ...]]
        ),
        vertex_alignments: dict[int, int] | None = None,
        vertex_positions: dict[int, np.ndarray] | None = None,
        reaction_factory: ReactionFactory = GenericReactionFactory(),
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
        scale_multiplier: float = 1.0,
    ) -> None:
        """
        Parameters:

            corners:
                Can be a :class:`dict` which maps the
                :class:`.BuildingBlock` instances to the ids of the
                vertices it should be placed on.

                Can also be a :class:`.BuildingBlock` instance, which
                should be placed on all corner vertices on the topology
                graph.

            linkers:
                Can be a :class:`dict` which maps the
                :class:`.BuildingBlock` instances to the ids of the
                vertices it should be placed on.

                Can also be a :class:`.BuildingBlock` instance, which
                should be placed on all linker vertices on the topology
                graph.

            vertex_alignments:
                A mapping from the id of a :class:`.Vertex`
                to an :class:`.Edge` connected to it.
                The :class:`.Edge` is used to align the first
                :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
                placed on that vertex. Only vertices which need to have
                their default edge changed need to be present in the
                :class:`dict`. If ``None`` then the default edge is used
                for each vertex. Changing which :class:`.Edge` is used will
                mean that the topology graph represents different
                structural isomers. The edge is referred to by a number
                between ``0`` (inclusive) and the number of edges the
                vertex is connected to (exclusive).

            vertex_positions:
                A mapping from the id of a :class:`.Vertex` to a custom
                :class:`.BuildingBlock` position. The default vertex
                alignment algorithm is still applied. Only vertices
                which need to have their default position changed need
                to be present in the :class:`dict`. Note that any
                vertices with modified positions will not be scaled like
                the rest of the building block positions and will not
                use neighbor placements in its positioning if requested
                by the default topology. If ``None`` then the default
                placement algorithm is used for each vertex.

            reaction_factory:
                The reaction factory to use for creating bonds between
                building blocks.

            num_processes:
                The number of parallel processes to create during
                :meth:`construct`.

            optimizer:
                Used to optimize the structure of the constructed
                molecule.

            scale_multiplier:
                Scales the positions of the vertices.

        Raises:

            :class:`AssertionError`
                If the any building block does not have a
                valid number of functional groups.

            :class:`ValueError`
                If the there are multiple building blocks with the
                same number of functional_groups in `building_blocks`,
                and they are not explicitly assigned to vertices. The
                desired placement of building blocks is ambiguous in
                this case.

            :class:`~.cage.UnoccupiedVertexError`
                If a vertex of the cage topology graph does not have a
                building block placed on it.

            :class:`~.cage.OverlyOccupiedVertexError`
                If a vertex of the cage topology graph has more than
                one building block placed on it.

        """

        if isinstance(corners, dict):
            building_blocks = corners
        else:
            building_blocks = {corners: (0, 1, 2, 3)}

        if isinstance(linkers, dict):
            linkers_dict = linkers
        else:
            linkers_dict = {linkers: (4, 5, 6, 7)}

        building_blocks.update(
            (building_block, vertices)
            for building_block, vertices in linkers_dict.items()
        )

        super().__init__(
            building_blocks,
            vertex_alignments=vertex_alignments,
            vertex_positions=vertex_positions,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
            scale_multiplier=scale_multiplier,
        )

    _vertex_prototypes = (
        AngledVertex(0, np.array([1, 1, 0])),
        AngledVertex(1, np.array([1, -1, 0])),
        AngledVertex(2, np.array([-1, -1, 0])),
        AngledVertex(3, np.array([-1, 1, 0])),
        LinearVertex(4, np.array([1, 0, 0]), False),
        LinearVertex(5, np.array([0, -1, 0]), False),
        LinearVertex(6, np.array([-1, 0, 0]), False),
        LinearVertex(7, np.array([0, 1, 0]), False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(2, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(4, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(6, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[7]),
    )

    _num_windows = 1
    _num_window_types = 1
