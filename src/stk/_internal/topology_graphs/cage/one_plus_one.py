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
from stk._internal.utilities.utilities import get_acute_vector

from .cage import Cage
from .vertices import NonLinearVertex


class OnePlusOneVertex(NonLinearVertex):
    def __init__(
        self,
        id,
        position,
        edge_normal,
        use_neighbor_placement=True,
        aligner_edge=0,
    ):
        super().__init__(
            id=id,
            position=position,
            use_neighbor_placement=use_neighbor_placement,
            aligner_edge=aligner_edge,
        )
        self._edge_normal = np.array(edge_normal, dtype=np.float64)

    def clone(self):
        clone = super().clone()
        clone._edge_normal = np.array(self._edge_normal)
        return clone

    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        edge_centroid = sum(edge.get_position() for edge in edges) / len(edges)
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        building_block = building_block.with_rotation_between_vectors(
            start=get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(
                    atom_ids=building_block.get_placer_ids(),
                ),
            ),
            target=self._edge_normal,
            origin=self._position,
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        start = fg_bonder_centroid - self._position
        edge_coord = edges[self._aligner_edge].get_position()
        building_block = building_block.with_rotation_to_minimize_angle(
            start=start,
            target=edge_coord - edge_centroid,
            axis=self._edge_normal,
            origin=self._position,
        )
        return building_block.get_position_matrix()


class OnePlusOne(Cage):
    """
    Represents a capsule cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='BrCC(CBr)CBr',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='Brc1cc(Br)cc(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.OnePlusOne({bb1: 0, bb2: 1}),
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    :class:`.Collapser` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='BrCC(CBr)CBr',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='Brc1cc(Br)cc(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.OnePlusOne(
                building_blocks={bb1: 0, bb2: 1},
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    Building blocks with three functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 1

    See :class:`.Cage` for more details and examples.

    """

    def __init__(
        self,
        building_blocks: (
            Iterable[BuildingBlock] | dict[BuildingBlock, tuple[int, ...]]
        ),
        vertex_alignments: dict[int, int] | None = None,
        reaction_factory: ReactionFactory = GenericReactionFactory(),
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
        scale_multiplier: float = 1.0,
    ) -> None:
        """
        Parameters:

            building_blocks:
                Can be a :class:`iterable` of :class:`.BuildingBlock`
                instances, which should be placed on the topology
                graph.

                Can also be a :class:`dict` which maps the
                :class:`.BuildingBlock` instances to the ids of the
                vertices it should be placed on. A :class:`dict` is
                required when there are multiple building blocks with
                the same number of functional groups, because in this
                case the desired placement is ambiguous.


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

        Notes:

            The `vertex_positions` argument is not available for this class
            due to strict positioning of the edges.

        """

        super().__init__(
            building_blocks,
            vertex_alignments=vertex_alignments,
            vertex_positions=None,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
            scale_multiplier=scale_multiplier,
        )

    _x = 1
    _vertex_prototypes = (
        OnePlusOneVertex(0, np.array([_x, 0.0, 0.0]), [1, 0, 0], False),
        OnePlusOneVertex(1, np.array([-_x, 0.0, 0.0]), [-1, 0, 0], False),
    )
    _edge_prototypes = (
        Edge(
            id=0,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            position=np.array([0.0, 1.0, 0.0]),
        ),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            position=np.array([0.0, -1.0, 1.0]),
        ),
        Edge(
            id=2,
            vertex1=_vertex_prototypes[0],
            vertex2=_vertex_prototypes[1],
            position=np.array([0.0, -1.0, -1.0]),
        ),
    )

    _num_windows = 3
    _num_window_types = 1
