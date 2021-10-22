"""
M3L3 Triangle
=============

"""

from __future__ import annotations

import numpy as np
import typing

from stk.utilities.typing import OneOrMany
from ..cage import Cage
from ..vertices import LinearVertex, AngledVertex
from ...topology_graph import NullOptimizer, Optimizer
from ....edge import Edge
from ....building_block import BuildingBlock
from ....reaction_factories import (
    ReactionFactory,
    GenericReactionFactory,
)


class M3L3Triangle(Cage):
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
                'C1=NC=CC(C2=CC=CC(C3=C'
                'C=NC=C3)=C2)=C1'
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
            topology_graph=stk.cage.M3L3Triangle(
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
                'C1=NC=CC(C2=CC=CC(C3=C'
                'C=NC=C3)=C2)=C1'
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
            topology_graph=stk.cage.M3L3Triangle(
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
        corners: typing.Union[
            BuildingBlock,
            dict[BuildingBlock, OneOrMany[int]]
        ],
        linkers: typing.Union[
            BuildingBlock,
            dict[BuildingBlock, OneOrMany[int]]
        ],
        vertex_alignments: typing.Optional[dict[int, int]] = None,
        reaction_factory: ReactionFactory = GenericReactionFactory(),
        num_processes: int = 1,
        optimizer: Optimizer = NullOptimizer(),
    ) -> None:
        """
        Initialize a :class:`.M3L3Triangle`.

        Parameters:

            corners:
                Can be a :class:`dict` which maps the
                :class:`.BuildingBlock` instances to the ids of the
                vertices it should be placed on. Valid numbers are
                0, 1 and 2.

                Can also be a :class:`.BuildingBlock` instance, which
                should be placed on all corner vertices on the topology
                graph.

            linkers:
                Can be a :class:`dict` which maps the
                :class:`.BuildingBlock` instances to the ids of the
                vertices it should be placed on. Valid numbers are
                3, 4 and 5.

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
                :class:`dict`. If ``None`` then the default edge is
                used for each vertex. Changing which :class:`.Edge` is
                used will mean that the topology graph represents
                different structural isomers. The edge is referred to
                by a number between ``0`` (inclusive) and the number of
                edges the vertex is connected to (exclusive).

            reaction_factory:
                The reaction factory to use for creating bonds between
                building blocks.

            num_processes:
                The number of parallel processes to create during
                :meth:`construct`.

            optimizer:
                Used to optimize the structure of the constructed
                molecule.

        """

        if isinstance(corners, dict):
            building_blocks = corners
        else:
            building_blocks = {corners: (0, 1, 2)}

        if isinstance(linkers, dict):
            linkers_dict = linkers
        else:
            linkers_dict = {linkers: (3, 4, 5)}

        building_blocks.update(
            (building_block, vertices)
            for building_block, vertices in linkers_dict.items()
        )

        super().__init__(
            building_blocks=building_blocks,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
        )

    _x = 2*np.sqrt(3)/4
    _y = 2
    _initial_vertex_prototypes = (
        AngledVertex(0, (0, _x, 0)),
        AngledVertex(1, (_y/2, -_x, 0)),
        AngledVertex(2, (-_y/2, -_x, 0)),
    )

    _vertex_prototypes = (
        *_initial_vertex_prototypes,

        LinearVertex.init_at_center(
            id=3,
            vertices=(
                _initial_vertex_prototypes[0],
                _initial_vertex_prototypes[1],
            ),
        ),
        LinearVertex.init_at_center(
            id=4,
            vertices=(
                _initial_vertex_prototypes[1],
                _initial_vertex_prototypes[2],
            ),
        ),
        LinearVertex.init_at_center(
            id=5,
            vertices=(
                _initial_vertex_prototypes[2],
                _initial_vertex_prototypes[0],
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[3]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[4]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[0], _vertex_prototypes[5]),
    )

    _num_windows = 1
    _num_window_types = 1

    def clone(self) -> M3L3Triangle:
        return self._clone()

    def with_building_blocks(
        self,
        building_block_map: dict[BuildingBlock, BuildingBlock]
    ) -> M3L3Triangle:
        return self._clone()._with_building_blocks(building_block_map)
