"""
M12L24
======

"""

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class M12L24(Cage):
    """
    Represents a cage topology graph.

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
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
            topology_graph=stk.cage.M12L24(
                building_blocks=(bb1, bb2),
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

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 11
        | 2-functional groups: 12 to 35

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [1, 0, 0]),
        NonLinearVertex(1, [-1, 0, 0]),
        NonLinearVertex(2, [0, 1, 0]),
        NonLinearVertex(3, [0, -1, 0]),
        NonLinearVertex(4, [0.5, 0.5, 0.707]),
        NonLinearVertex(5, [0.5, -0.5, 0.707]),
        NonLinearVertex(6, [-0.5, 0.5, 0.707]),
        NonLinearVertex(7, [-0.5, -0.5, 0.707]),
        NonLinearVertex(8, [0.5, 0.5, -0.707]),
        NonLinearVertex(9, [0.5, -0.5, -0.707]),
        NonLinearVertex(10, [-0.5, 0.5, -0.707]),
        NonLinearVertex(11, [-0.5, -0.5, -0.707]),

        LinearVertex(12, [0.9, 0.31, 0.31], False),
        LinearVertex(13, [0.9, 0.31, -0.31], False),
        LinearVertex(14, [0.9, -0.31, 0.31], False),
        LinearVertex(15, [0.9, -0.31, -0.31], False),

        LinearVertex(16, [-0.9, 0.31, 0.31], False),
        LinearVertex(17, [-0.9, 0.31, -0.31], False),
        LinearVertex(18, [-0.9, -0.31, 0.31], False),
        LinearVertex(19, [-0.9, -0.31, -0.31], False),

        LinearVertex(20, [0.31, 0.9, 0.31], False),
        LinearVertex(21, [0.31, 0.9, -0.31], False),
        LinearVertex(22, [-0.31, 0.9, 0.31], False),
        LinearVertex(23, [-0.31, 0.9, -0.31], False),

        LinearVertex(24, [0.31, -0.9, 0.31], False),
        LinearVertex(25, [0.31, -0.9, -0.31], False),
        LinearVertex(26, [-0.31, -0.9, 0.31], False),
        LinearVertex(27, [-0.31, -0.9, -0.31], False),

        LinearVertex(28, [0.58, 0, 0.82], False),
        LinearVertex(29, [-0.58, 0, 0.82], False),
        LinearVertex(30, [0, 0.58, 0.82], False),
        LinearVertex(31, [0, -0.58, 0.82], False),
        LinearVertex(32, [0.58, 0, -0.82], False),
        LinearVertex(33, [-0.58, 0, -0.82], False),
        LinearVertex(34, [0, 0.58, -0.82], False),
        LinearVertex(35, [0, -0.58, -0.82], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[12]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[13]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[14]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[15]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[16]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[17]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[18]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[19]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[20]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[21]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[22]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[23]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[24]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[25]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[26]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[27]),

        Edge(16, _vertex_prototypes[4], _vertex_prototypes[28]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[30]),
        Edge(18, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(19, _vertex_prototypes[4], _vertex_prototypes[20]),

        Edge(20, _vertex_prototypes[5], _vertex_prototypes[14]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[24]),
        Edge(22, _vertex_prototypes[5], _vertex_prototypes[28]),
        Edge(23, _vertex_prototypes[5], _vertex_prototypes[31]),

        Edge(24, _vertex_prototypes[6], _vertex_prototypes[16]),
        Edge(25, _vertex_prototypes[6], _vertex_prototypes[29]),
        Edge(26, _vertex_prototypes[6], _vertex_prototypes[30]),
        Edge(27, _vertex_prototypes[6], _vertex_prototypes[22]),

        Edge(28, _vertex_prototypes[7], _vertex_prototypes[18]),
        Edge(29, _vertex_prototypes[7], _vertex_prototypes[26]),
        Edge(30, _vertex_prototypes[7], _vertex_prototypes[31]),
        Edge(31, _vertex_prototypes[7], _vertex_prototypes[29]),

        Edge(32, _vertex_prototypes[8], _vertex_prototypes[13]),
        Edge(33, _vertex_prototypes[8], _vertex_prototypes[32]),
        Edge(34, _vertex_prototypes[8], _vertex_prototypes[34]),
        Edge(35, _vertex_prototypes[8], _vertex_prototypes[21]),

        Edge(36, _vertex_prototypes[9], _vertex_prototypes[15]),
        Edge(37, _vertex_prototypes[9], _vertex_prototypes[32]),
        Edge(38, _vertex_prototypes[9], _vertex_prototypes[35]),
        Edge(39, _vertex_prototypes[9], _vertex_prototypes[25]),

        Edge(40, _vertex_prototypes[10], _vertex_prototypes[17]),
        Edge(41, _vertex_prototypes[10], _vertex_prototypes[23]),
        Edge(42, _vertex_prototypes[10], _vertex_prototypes[34]),
        Edge(43, _vertex_prototypes[10], _vertex_prototypes[33]),

        Edge(44, _vertex_prototypes[11], _vertex_prototypes[19]),
        Edge(45, _vertex_prototypes[11], _vertex_prototypes[33]),
        Edge(46, _vertex_prototypes[11], _vertex_prototypes[27]),
        Edge(47, _vertex_prototypes[11], _vertex_prototypes[35]),
    )

    _num_windows = 14
    _num_window_types = 2
