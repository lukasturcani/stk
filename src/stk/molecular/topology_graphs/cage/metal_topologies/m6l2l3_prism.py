"""
M6L2L3 Prism
============

"""

import numpy as np

from ..cage import Cage
from ..vertices import NonLinearVertex
from ...topology_graph import Edge


class M6L2L3Prism(Cage):
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
                for i in range(3)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles='c2cnccc2c1cc(c4ccncc4)cc(c3ccncc3)c1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        bb3 = stk.BuildingBlock(
            smiles=(
                'c2cnccc2C#Cc1c(C#Cc4ccncc4)cc'
                '(C#Cc3ccncc3)c(C#Cc5ccncc5)c1'
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
            topology_graph=stk.cage.M6L2L3Prism(
                building_blocks={
                    bb1: range(6),
                    bb2: (6, 7),
                    bb3: range(8, 11),
                },
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

    :class:`.Collapser` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(3)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
            smiles='c2cnccc2c1cc(c4ccncc4)cc(c3ccncc3)c1',
            functional_groups=[
                stk.SmartsFunctionalGroupFactory(
                    smarts='[#6]~[#7X2]~[#6]',
                    bonders=(1, ),
                    deleters=(),
                ),
            ],
        )

        bb3 = stk.BuildingBlock(
            smiles=(
                'c2cnccc2C#Cc1c(C#Cc4ccncc4)cc'
                '(C#Cc3ccncc3)c(C#Cc5ccncc5)c1'
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
            topology_graph=stk.cage.M6L2L3Prism(
                building_blocks={
                    bb1: range(6),
                    bb2: (6, 7),
                    bb3: range(8, 11),
                },
                optimizer=stk.Collapser(),
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

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with three and four functional groups are
    required for ligand type A and B, respectively, on this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups (metal): 0 to 5
        | 3-functional groups (ligand A): 6, 7
        | 4-functional groups (ligand B): 8 to 10

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [-1, -1/np.sqrt(3), 1]),
        NonLinearVertex(1, [1, -1/np.sqrt(3), 1]),
        NonLinearVertex(2, [0, 2/np.sqrt(3), 1]),

        NonLinearVertex(3, [-1, -1/np.sqrt(3), -1]),
        NonLinearVertex(4, [1, -1/np.sqrt(3), -1]),
        NonLinearVertex(5, [0, 2/np.sqrt(3), -1]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        NonLinearVertex.init_at_center(
            id=6,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[2],
            ),
        ),
        NonLinearVertex.init_at_center(
            id=7,
            vertices=(
                _vertex_prototypes[3],
                _vertex_prototypes[4],
                _vertex_prototypes[5],
            ),
        ),

        NonLinearVertex.init_at_center(
            id=8,
            vertices=(
                _vertex_prototypes[0],
                _vertex_prototypes[1],
                _vertex_prototypes[3],
                _vertex_prototypes[4],
            ),
        ),
        NonLinearVertex.init_at_center(
            id=9,
            vertices=(
                _vertex_prototypes[1],
                _vertex_prototypes[2],
                _vertex_prototypes[4],
                _vertex_prototypes[5],
            ),
        ),
        NonLinearVertex.init_at_center(
            id=10,
            vertices=(
                _vertex_prototypes[2],
                _vertex_prototypes[0],
                _vertex_prototypes[5],
                _vertex_prototypes[3],
            ),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[10]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[9]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(7, _vertex_prototypes[2], _vertex_prototypes[9]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[3], _vertex_prototypes[10]),
        Edge(12, _vertex_prototypes[4], _vertex_prototypes[7]),
        Edge(13, _vertex_prototypes[4], _vertex_prototypes[8]),
        Edge(14, _vertex_prototypes[4], _vertex_prototypes[9]),
        Edge(15, _vertex_prototypes[5], _vertex_prototypes[7]),
        Edge(16, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(17, _vertex_prototypes[5], _vertex_prototypes[10]),
    )

    _num_windows = 4
    _num_window_types = 1
