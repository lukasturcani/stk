"""
M6L8L12 Cuboctahedron
==========

"""

import numpy as np

from stk._internal.topology_graphs.edge import Edge

from .cage import Cage
from .vertices import LinearVertex, NonLinearVertex


class M6L6L12Cuboctahedron(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        m1 = stk.BuildingBlock('[Ti](Br)(Br)(Br)(Br)(Br)(Br)(Br)(Br)', functional_groups=[stk.BromoFactory()])

        bb1 = stk.BuildingBlock('C(OBr)=[O+]Br', functional_groups=[stk.BromoFactory()])

        bb2 = stk.BuildingBlock('[O+](Br)(Br)(Br)', functional_groups=[stk.BromoFactory()])

        cage = stk.ConstructedMolecule(
            topology_graph = stk.cage.M6L6L12Cuboctahedron(
                building_blocks = {
                    m1:range(0,6),
                    bb1:range(6,18),
                    bb2:range(18,26)
                },
                optimizer=stk.Nulloptimizer(),
            )
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

        m1 = stk.BuildingBlock('[Ti](Br)(Br)(Br)(Br)(Br)(Br)(Br)(Br)', functional_groups=[stk.BromoFactory()])

        bb1 = stk.BuildingBlock('C(OBr)=[O+]Br', functional_groups=[stk.BromoFactory()])

        bb2 = stk.BuildingBlock('[O+](Br)(Br)(Br)', functional_groups=[stk.BromoFactory()])

        cage = stk.ConstructedMolecule(
            topology_graph = stk.cage.M6L6L12Cuboctahedron(
                building_blocks = {
                    m1:range(0,6),
                    bb1:range(6,18),
                    bb2:range(18,26)
                },
                optimizer=stk.Collapser(),
            )
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

    Metal building blocks with eight functional groups are
    required for this topology.

    One ligand building blocks with two functional groups and another ligand 
    building block with three functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 8-functional groups: 0  to 5
        | 2-functional groups: 6  to 17
        | 3-functional groups: 18 to 25

    See :class:`.Cage` for more details and examples.

    """

    _x = np.sqrt(2)
    _vertex_prototypes = (
        # M
        NonLinearVertex(0, np.array([1, 0, 0])),
        NonLinearVertex(1, np.array([0, 1, 0])),
        NonLinearVertex(2, np.array([-1, 0, 0])),
        NonLinearVertex(3, np.array([0, -1, 0])),
        NonLinearVertex(4, np.array([0, 0, 1])),
        NonLinearVertex(5, np.array([0, 0, -1])),
        # L1
        LinearVertex(6, np.array([2, 2, 0]), False),
        LinearVertex(7, np.array([2, -2, 0]), False),
        LinearVertex(8, np.array([2, 0, 2]), False),
        LinearVertex(9, np.array([2, 0, -2]), False),
        LinearVertex(10, np.array([-2, 2, 0]), False),
        LinearVertex(11, np.array([-2, -2, 0]), False),
        LinearVertex(12, np.array([-2, 0, 2]), False),
        LinearVertex(13, np.array([-2, 0, -2]), False),
        LinearVertex(14, np.array([0, 2, 2]), False),
        LinearVertex(15, np.array([0, 2, -2]), False),
        LinearVertex(16, np.array([0, -2, 2]), False),
        LinearVertex(17, np.array([0, -2, -2]), False),
        # L2
        LinearVertex(18, np.array([ _x, _x, _x]), False),
        LinearVertex(19, np.array([-_x, _x, _x]), False),
        LinearVertex(20, np.array([-_x,-_x, _x]), False),
        LinearVertex(21, np.array([ _x,-_x, _x]), False),
        LinearVertex(22, np.array([ _x, _x,-_x]), False),
        LinearVertex(23, np.array([-_x, _x,-_x]), False),
        LinearVertex(24, np.array([-_x,-_x,-_x]), False),
        LinearVertex(25, np.array([ _x,-_x,-_x]), False),
    )

    _edge_prototypes = (
        # Outer shell
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[7]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[9]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[10]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[14]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[15]),
        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[12]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[13]),
        Edge(12, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[11]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[16]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[17]),
        Edge(16, _vertex_prototypes[4], _vertex_prototypes[8]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(18, _vertex_prototypes[4], _vertex_prototypes[14]),
        Edge(19, _vertex_prototypes[4], _vertex_prototypes[16]),
        Edge(20, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[13]),
        Edge(22, _vertex_prototypes[5], _vertex_prototypes[15]),
        Edge(23, _vertex_prototypes[5], _vertex_prototypes[17]),
        # Inner shell
        Edge(24, _vertex_prototypes[0], _vertex_prototypes[18]),
        Edge(25, _vertex_prototypes[0], _vertex_prototypes[21]),
        Edge(26, _vertex_prototypes[0], _vertex_prototypes[22]),
        Edge(27, _vertex_prototypes[0], _vertex_prototypes[25]),
        Edge(28, _vertex_prototypes[1], _vertex_prototypes[18]),
        Edge(29, _vertex_prototypes[1], _vertex_prototypes[19]),
        Edge(30, _vertex_prototypes[1], _vertex_prototypes[22]),
        Edge(31, _vertex_prototypes[1], _vertex_prototypes[23]),
        Edge(32, _vertex_prototypes[2], _vertex_prototypes[19]),
        Edge(33, _vertex_prototypes[2], _vertex_prototypes[20]),
        Edge(34, _vertex_prototypes[2], _vertex_prototypes[23]),
        Edge(35, _vertex_prototypes[2], _vertex_prototypes[24]),
        Edge(36, _vertex_prototypes[3], _vertex_prototypes[20]),
        Edge(37, _vertex_prototypes[3], _vertex_prototypes[21]),
        Edge(38, _vertex_prototypes[3], _vertex_prototypes[24]),
        Edge(39, _vertex_prototypes[3], _vertex_prototypes[25]),
        Edge(40, _vertex_prototypes[4], _vertex_prototypes[18]),
        Edge(41, _vertex_prototypes[4], _vertex_prototypes[19]),
        Edge(42, _vertex_prototypes[4], _vertex_prototypes[20]),
        Edge(43, _vertex_prototypes[4], _vertex_prototypes[21]),
        Edge(44, _vertex_prototypes[5], _vertex_prototypes[22]),
        Edge(45, _vertex_prototypes[5], _vertex_prototypes[23]),
        Edge(46, _vertex_prototypes[5], _vertex_prototypes[24]),
        Edge(47, _vertex_prototypes[5], _vertex_prototypes[25]),
    )

    _num_windows = 8
    _num_window_types = 1
