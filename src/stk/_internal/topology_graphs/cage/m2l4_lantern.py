"""
M2L4 Lantern
============

"""

from stk._internal.topology_graphs.edge import Edge

from .cage import Cage
from .vertices import LinearVertex, NonLinearVertex


class M2L4Lantern(Cage):
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
            topology_graph=stk.cage.M2L4Lantern(
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

    :class:`.MCHammer` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        four_c_bb = stk.BuildingBlock(
            smiles="[Br][C]([Br])([Br])[Br]",
            position_matrix=[
                [-2, 0, -1],
                [0, 0, 1],
                [0, -2, -1],
                [2, 0, 1],
                [0, 2, 1],
            ],
            functional_groups=(stk.BromoFactory(placers=(0, 1)),),
        )

        two_c_bb = stk.BuildingBlock(
            smiles="[Br][N][Br]",
            position_matrix=[
                [-2, 0, -1],
                [0, 0, 1],
                [0, -2, -1],
            ],
            functional_groups=(stk.BromoFactory(placers=(0, 1)),),
        )

        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.M2L4Lantern(
                building_blocks=(four_c_bb, two_c_bb),
            ),
        )
        rgb1 = [192, 87, 161]
        rgb2 = [97, 201, 217]
        points = cage.get_num_atoms()
        colour_list = [
            rgb1[0] + ((rgb2[0]-rgb1[0])/points)*i,
            rgb1[1] + ((rgb2[1]-rgb1[1])/points)*i,
            rgb1[2] + ((rgb2[2]-rgb1[2])/points)*i,
            for i in range(points)
        ]

        moldoc_display_molecule = molecule.Molecule(
            atoms=(
                molecule.Atom(
                    atomic_number=atom.get_atomic_number(),
                    position=position,
                    config=molecule.AtomConfig(
                        color=molecule.Color(
                            red=col[0],
                            green=col[1],
                            blue=col[2],
                        ),
                        size=2.0,
                    )
                ) for atom, position, col in zip(
                    cage.get_atoms(),
                    cage.get_position_matrix(),
                    colour_list,
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

    Metal building blocks (COLOUR1) with four functional groups
    are required for this topology.

    Ligand building blocks (COLOUR2) with two functional groups
    are required for this topology.

    :class:`.MCHammer` optimization is recommended for construction.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups (COLOUR1): 0 to 1
        | 2-functional groups (COLOUR2): 2 to 5

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [0, 0.5, 0]),
        NonLinearVertex(1, [0, -0.5, 0]),
        LinearVertex(2, [1, 0, 0], False),
        LinearVertex(3, [0, 0, 1], False),
        LinearVertex(4, [-1, 0, 0], False),
        LinearVertex(5, [0, 0, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[5]),
    )

    _num_windows = 4
    _num_window_types = 1
