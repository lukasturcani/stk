"""
Ten Plus Twenty
===============

"""

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class TenPlusTwenty(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='BrCCBr',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='Brc1c(Br)cc(Br)c(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.TenPlusTwenty((bb1, bb2)),
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
            smiles='BrCCBr',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='Brc1c(Br)cc(Br)c(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.TenPlusTwenty(
                building_blocks=(bb1, bb2),
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
                    order=bond.get_order(),
                ) for bond in cage.get_bonds()
            ),
        )

    Nonlinear building blocks with four functional groups are
    required for this topology.

    Linear building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 9
        | 2-functional groups: 10 to 29

    See :class:`.Cage` for more details and examples.

    """

    _x = 1.5
    _z = _x / 2
    _non_linears = (
        NonLinearVertex(0, [-_x, _x, -_z]),
        NonLinearVertex(1, [-_x, -_x, -_z]),
        NonLinearVertex(2, [_x, _x, -_z]),
        NonLinearVertex(3, [_x, -_x, -_z]),

        NonLinearVertex(4, [-_x, _x, _z]),
        NonLinearVertex(5, [-_x, -_x, _z]),
        NonLinearVertex(6, [_x, _x, _z]),
        NonLinearVertex(7, [_x, -_x, _z]),

        NonLinearVertex(8, [0, 0, _x]),
        NonLinearVertex(9, [0, 0, -_x]),
    )

    _vertex_prototypes = (
        *_non_linears,

        LinearVertex.init_at_center(
            id=10,
            vertices=(_non_linears[0], _non_linears[2]),
        ),
        LinearVertex.init_at_center(
            id=11,
            vertices=(_non_linears[0], _non_linears[1]),
        ),
        LinearVertex.init_at_center(
            id=12,
            vertices=(_non_linears[1], _non_linears[3]),
        ),
        LinearVertex.init_at_center(
            id=13,
            vertices=(_non_linears[2], _non_linears[3]),
        ),

        LinearVertex.init_at_center(
            id=14,
            vertices=(_non_linears[4], _non_linears[6]),
        ),
        LinearVertex.init_at_center(
            id=15,
            vertices=(_non_linears[4], _non_linears[5]),
        ),
        LinearVertex.init_at_center(
            id=16,
            vertices=(_non_linears[5], _non_linears[7]),
        ),
        LinearVertex.init_at_center(
            id=17,
            vertices=(_non_linears[6], _non_linears[7]),
        ),

        LinearVertex.init_at_center(
            id=18,
            vertices=(_non_linears[0], _non_linears[4]),
        ),
        LinearVertex.init_at_center(
            id=19,
            vertices=(_non_linears[1], _non_linears[5]),
        ),
        LinearVertex.init_at_center(
            id=20,
            vertices=(_non_linears[2], _non_linears[6]),
        ),
        LinearVertex.init_at_center(
            id=21,
            vertices=(_non_linears[3], _non_linears[7]),
        ),

        LinearVertex.init_at_center(
            id=22,
            vertices=(_non_linears[8], _non_linears[4]),
        ),
        LinearVertex.init_at_center(
            id=23,
            vertices=(_non_linears[8], _non_linears[5]),
        ),
        LinearVertex.init_at_center(
            id=24,
            vertices=(_non_linears[8], _non_linears[6]),
        ),
        LinearVertex.init_at_center(
            id=25,
            vertices=(_non_linears[8], _non_linears[7]),
        ),

        LinearVertex.init_at_center(
            id=26,
            vertices=(_non_linears[9], _non_linears[0]),
        ),
        LinearVertex.init_at_center(
            id=27,
            vertices=(_non_linears[9], _non_linears[1]),
        ),
        LinearVertex.init_at_center(
            id=28,
            vertices=(_non_linears[9], _non_linears[2]),
        ),
        LinearVertex.init_at_center(
            id=29,
            vertices=(_non_linears[9], _non_linears[3]),
        )
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[10], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[10], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[11], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[11], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[12], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[12], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[13], _vertex_prototypes[2]),
        Edge(7, _vertex_prototypes[13], _vertex_prototypes[3]),

        Edge(8, _vertex_prototypes[14], _vertex_prototypes[4]),
        Edge(9, _vertex_prototypes[14], _vertex_prototypes[6]),
        Edge(10, _vertex_prototypes[15], _vertex_prototypes[4]),
        Edge(11, _vertex_prototypes[15], _vertex_prototypes[5]),

        Edge(12, _vertex_prototypes[16], _vertex_prototypes[5]),
        Edge(13, _vertex_prototypes[16], _vertex_prototypes[7]),
        Edge(14, _vertex_prototypes[17], _vertex_prototypes[6]),
        Edge(15, _vertex_prototypes[17], _vertex_prototypes[7]),

        Edge(16, _vertex_prototypes[18], _vertex_prototypes[0]),
        Edge(17, _vertex_prototypes[18], _vertex_prototypes[4]),
        Edge(18, _vertex_prototypes[19], _vertex_prototypes[1]),
        Edge(19, _vertex_prototypes[19], _vertex_prototypes[5]),

        Edge(20, _vertex_prototypes[20], _vertex_prototypes[2]),
        Edge(21, _vertex_prototypes[20], _vertex_prototypes[6]),
        Edge(22, _vertex_prototypes[21], _vertex_prototypes[3]),
        Edge(23, _vertex_prototypes[21], _vertex_prototypes[7]),

        Edge(24, _vertex_prototypes[22], _vertex_prototypes[8]),
        Edge(25, _vertex_prototypes[22], _vertex_prototypes[4]),
        Edge(26, _vertex_prototypes[23], _vertex_prototypes[8]),
        Edge(27, _vertex_prototypes[23], _vertex_prototypes[5]),

        Edge(28, _vertex_prototypes[24], _vertex_prototypes[8]),
        Edge(29, _vertex_prototypes[24], _vertex_prototypes[6]),
        Edge(30, _vertex_prototypes[25], _vertex_prototypes[8]),
        Edge(31, _vertex_prototypes[25], _vertex_prototypes[7]),

        Edge(32, _vertex_prototypes[26], _vertex_prototypes[9]),
        Edge(33, _vertex_prototypes[26], _vertex_prototypes[0]),
        Edge(34, _vertex_prototypes[27], _vertex_prototypes[9]),
        Edge(35, _vertex_prototypes[27], _vertex_prototypes[1]),

        Edge(36, _vertex_prototypes[28], _vertex_prototypes[9]),
        Edge(37, _vertex_prototypes[28], _vertex_prototypes[2]),
        Edge(38, _vertex_prototypes[29], _vertex_prototypes[9]),
        Edge(39, _vertex_prototypes[29], _vertex_prototypes[3]),
    )

    _num_windows = 12
    _num_window_types = 2
