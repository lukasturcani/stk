"""
Twenty Plus Thirty
==================

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class TwentyPlusThirty(Cage):
    """
    Represents a cage topology graph.

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='BrCCBr',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='Brc1cc(Br)cc(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.TwentyPlusThirty((bb1, bb2)),
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

    Nonlinear building blocks with three functional groups are
    required for this topology.

    Linear building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 19
        | 2-functional groups: 20 to 49

    See :class:`.Cage` for more details and examples.

    """

    # Source: http://tinyurl.com/h2dl949
    _phi = (1 + np.sqrt(5))/2
    _x = 1.5
    _non_linears = (
        NonLinearVertex(0, [_x*_phi, 0.0, _x/_phi]),
        NonLinearVertex(1, [_x*-_phi, 0.0, _x/_phi]),
        NonLinearVertex(2, [_x*-_phi, 0.0, _x/-_phi]),
        NonLinearVertex(3, [_x*_phi, 0.0, _x/-_phi]),

        NonLinearVertex(4, [_x/_phi, _x*_phi, 0.0]),
        NonLinearVertex(5, [_x/_phi, _x*-_phi, 0.0]),
        NonLinearVertex(6, [_x/-_phi, _x*-_phi, 0.0]),
        NonLinearVertex(7, [_x/-_phi, _x*_phi, 0.0]),
        NonLinearVertex(8, [0.0, _x/_phi, _x*_phi]),
        NonLinearVertex(9, [0.0, _x/_phi, _x*-_phi]),
        NonLinearVertex(10, [0.0, _x/-_phi, _x*-_phi]),
        NonLinearVertex(11, [0.0, _x/-_phi, _x*_phi]),

        NonLinearVertex(12, [_x, _x, _x]),
        NonLinearVertex(13, [_x, -_x, _x]),
        NonLinearVertex(14, [-_x, -_x, _x]),
        NonLinearVertex(15, [-_x, _x, _x]),
        NonLinearVertex(16, [-_x, _x, -_x]),
        NonLinearVertex(17, [_x, _x, -_x]),
        NonLinearVertex(18, [_x, -_x, -_x]),
        NonLinearVertex(19, [-_x, -_x, -_x]),
    )

    _vertex_prototypes = (
        *_non_linears,
        LinearVertex.init_at_center(
            id=20,
            vertices=(_non_linears[0], _non_linears[13]),
        ),
        LinearVertex.init_at_center(
            id=21,
            vertices=(_non_linears[0], _non_linears[12]),
        ),
        LinearVertex.init_at_center(
            id=22,
            vertices=(_non_linears[0], _non_linears[3]),
        ),

        LinearVertex.init_at_center(
            id=23,
            vertices=(_non_linears[1], _non_linears[14]),
        ),
        LinearVertex.init_at_center(
            id=24,
            vertices=(_non_linears[1], _non_linears[15]),
        ),
        LinearVertex.init_at_center(
            id=25,
            vertices=(_non_linears[1], _non_linears[2]),
        ),

        LinearVertex.init_at_center(
            id=26,
            vertices=(_non_linears[2], _non_linears[19]),
        ),
        LinearVertex.init_at_center(
            id=27,
            vertices=(_non_linears[2], _non_linears[16]),
        ),

        LinearVertex.init_at_center(
            id=28,
            vertices=(_non_linears[3], _non_linears[18]),
        ),
        LinearVertex.init_at_center(
            id=29,
            vertices=(_non_linears[3], _non_linears[17]),
        ),

        LinearVertex.init_at_center(
            id=30,
            vertices=(_non_linears[4], _non_linears[12]),
        ),
        LinearVertex.init_at_center(
            id=31,
            vertices=(_non_linears[4], _non_linears[7]),
        ),
        LinearVertex.init_at_center(
            id=32,
            vertices=(_non_linears[4], _non_linears[17]),
        ),

        LinearVertex.init_at_center(
            id=33,
            vertices=(_non_linears[5], _non_linears[6]),
        ),
        LinearVertex.init_at_center(
            id=34,
            vertices=(_non_linears[5], _non_linears[18]),
        ),
        LinearVertex.init_at_center(
            id=35,
            vertices=(_non_linears[5], _non_linears[13]),
        ),

        LinearVertex.init_at_center(
            id=36,
            vertices=(_non_linears[6], _non_linears[14]),
        ),
        LinearVertex.init_at_center(
            id=37,
            vertices=(_non_linears[6], _non_linears[19]),
        ),

        LinearVertex.init_at_center(
            id=38,
            vertices=(_non_linears[7], _non_linears[15]),
        ),
        LinearVertex.init_at_center(
            id=39,
            vertices=(_non_linears[7], _non_linears[16]),
        ),

        LinearVertex.init_at_center(
            id=40,
            vertices=(_non_linears[8], _non_linears[11]),
        ),
        LinearVertex.init_at_center(
            id=41,
            vertices=(_non_linears[8], _non_linears[12]),
        ),
        LinearVertex.init_at_center(
            id=42,
            vertices=(_non_linears[8], _non_linears[15]),
        ),

        LinearVertex.init_at_center(
            id=43,
            vertices=(_non_linears[9], _non_linears[10]),
        ),
        LinearVertex.init_at_center(
            id=44,
            vertices=(_non_linears[9], _non_linears[17]),
        ),
        LinearVertex.init_at_center(
            id=45,
            vertices=(_non_linears[9], _non_linears[16]),
        ),

        LinearVertex.init_at_center(
            id=46,
            vertices=(_non_linears[10], _non_linears[18]),
        ),
        LinearVertex.init_at_center(
            id=47,
            vertices=(_non_linears[10], _non_linears[19]),
        ),

        LinearVertex.init_at_center(
            id=48,
            vertices=(_non_linears[11], _non_linears[14]),
        ),
        LinearVertex.init_at_center(
            id=49,
            vertices=(_non_linears[11], _non_linears[13])),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[20], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[20], _vertex_prototypes[13]),
        Edge(2, _vertex_prototypes[21], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[21], _vertex_prototypes[12]),
        Edge(4, _vertex_prototypes[22], _vertex_prototypes[0]),
        Edge(5, _vertex_prototypes[22], _vertex_prototypes[3]),

        Edge(6, _vertex_prototypes[23], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[23], _vertex_prototypes[14]),
        Edge(8, _vertex_prototypes[24], _vertex_prototypes[1]),
        Edge(9, _vertex_prototypes[24], _vertex_prototypes[15]),
        Edge(10, _vertex_prototypes[25], _vertex_prototypes[1]),
        Edge(11, _vertex_prototypes[25], _vertex_prototypes[2]),

        Edge(12, _vertex_prototypes[26], _vertex_prototypes[2]),
        Edge(13, _vertex_prototypes[26], _vertex_prototypes[19]),
        Edge(14, _vertex_prototypes[27], _vertex_prototypes[2]),
        Edge(15, _vertex_prototypes[27], _vertex_prototypes[16]),

        Edge(16, _vertex_prototypes[28], _vertex_prototypes[3]),
        Edge(17, _vertex_prototypes[28], _vertex_prototypes[18]),
        Edge(18, _vertex_prototypes[29], _vertex_prototypes[3]),
        Edge(19, _vertex_prototypes[29], _vertex_prototypes[17]),

        Edge(20, _vertex_prototypes[30], _vertex_prototypes[4]),
        Edge(21, _vertex_prototypes[30], _vertex_prototypes[12]),
        Edge(22, _vertex_prototypes[31], _vertex_prototypes[4]),
        Edge(23, _vertex_prototypes[31], _vertex_prototypes[7]),
        Edge(24, _vertex_prototypes[32], _vertex_prototypes[4]),
        Edge(25, _vertex_prototypes[32], _vertex_prototypes[17]),

        Edge(26, _vertex_prototypes[33], _vertex_prototypes[5]),
        Edge(27, _vertex_prototypes[33], _vertex_prototypes[6]),
        Edge(28, _vertex_prototypes[34], _vertex_prototypes[5]),
        Edge(29, _vertex_prototypes[34], _vertex_prototypes[18]),
        Edge(30, _vertex_prototypes[35], _vertex_prototypes[5]),
        Edge(31, _vertex_prototypes[35], _vertex_prototypes[13]),

        Edge(32, _vertex_prototypes[36], _vertex_prototypes[6]),
        Edge(33, _vertex_prototypes[36], _vertex_prototypes[14]),
        Edge(34, _vertex_prototypes[37], _vertex_prototypes[6]),
        Edge(35, _vertex_prototypes[37], _vertex_prototypes[19]),

        Edge(36, _vertex_prototypes[38], _vertex_prototypes[7]),
        Edge(37, _vertex_prototypes[38], _vertex_prototypes[15]),
        Edge(38, _vertex_prototypes[39], _vertex_prototypes[7]),
        Edge(39, _vertex_prototypes[39], _vertex_prototypes[16]),

        Edge(40, _vertex_prototypes[40], _vertex_prototypes[8]),
        Edge(41, _vertex_prototypes[40], _vertex_prototypes[11]),
        Edge(42, _vertex_prototypes[41], _vertex_prototypes[8]),
        Edge(43, _vertex_prototypes[41], _vertex_prototypes[12]),
        Edge(44, _vertex_prototypes[42], _vertex_prototypes[8]),
        Edge(45, _vertex_prototypes[42], _vertex_prototypes[15]),

        Edge(46, _vertex_prototypes[43], _vertex_prototypes[9]),
        Edge(47, _vertex_prototypes[43], _vertex_prototypes[10]),
        Edge(48, _vertex_prototypes[44], _vertex_prototypes[9]),
        Edge(49, _vertex_prototypes[44], _vertex_prototypes[17]),
        Edge(50, _vertex_prototypes[45], _vertex_prototypes[9]),
        Edge(51, _vertex_prototypes[45], _vertex_prototypes[16]),

        Edge(52, _vertex_prototypes[46], _vertex_prototypes[10]),
        Edge(53, _vertex_prototypes[46], _vertex_prototypes[18]),
        Edge(54, _vertex_prototypes[47], _vertex_prototypes[10]),
        Edge(55, _vertex_prototypes[47], _vertex_prototypes[19]),

        Edge(56, _vertex_prototypes[48], _vertex_prototypes[11]),
        Edge(57, _vertex_prototypes[48], _vertex_prototypes[14]),
        Edge(58, _vertex_prototypes[49], _vertex_prototypes[11]),
        Edge(59, _vertex_prototypes[49], _vertex_prototypes[13]),
    )

    _num_windows = 12
    _num_window_types = 1
