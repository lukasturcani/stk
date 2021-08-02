"""
M8L6 Cube
=========

"""

from ..cage import Cage
from ..vertices import NonLinearVertex
from ...topology_graph import Edge


class M8L6Cube(Cage):
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
            topology_graph=stk.cage.M8L6Cube(
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

        bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(3)
            ),
            position_matrix=[[0, 0, 0]],
        )

        bb2 = stk.BuildingBlock(
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
            topology_graph=stk.cage.M8L6Cube(
                building_blocks=(bb1, bb2),
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

    Metal building blocks with three functional groups are
    required for this topology.

    Ligand building blocks with four functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 7
        | 4-functional groups: 8 to 13

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(
            id=0,
            position=[1, 1, 1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=1,
            position=[1, -1, 1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=2,
            position=[-1, -1, 1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=3,
            position=[-1, 1, 1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=4,
            position=[1, 1, -1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=5,
            position=[1, -1, -1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=6,
            position=[-1, -1, -1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=7,
            position=[-1, 1, -1],
            use_neighbor_placement=False,
        ),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        NonLinearVertex(
            id=8,
            position=[0, 0, 1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=9,
            position=[1, 0, 0],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=10,
            position=[0, 1, 0],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=11,
            position=[-1, 0, 0],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=12,
            position=[0, 0, -1],
            use_neighbor_placement=False,
        ),
        NonLinearVertex(
            id=13,
            position=[0, -1, 0],
            use_neighbor_placement=False,
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[8]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[8]),
        Edge(2, _vertex_prototypes[2], _vertex_prototypes[8]),
        Edge(3, _vertex_prototypes[3], _vertex_prototypes[8]),

        Edge(4, _vertex_prototypes[4], _vertex_prototypes[9]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[9]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[9]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[9]),

        Edge(8, _vertex_prototypes[4], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[0], _vertex_prototypes[10]),
        Edge(10, _vertex_prototypes[3], _vertex_prototypes[10]),
        Edge(11, _vertex_prototypes[7], _vertex_prototypes[10]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[11]),
        Edge(13, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(14, _vertex_prototypes[6], _vertex_prototypes[11]),
        Edge(15, _vertex_prototypes[7], _vertex_prototypes[11]),

        Edge(16, _vertex_prototypes[5], _vertex_prototypes[12]),
        Edge(17, _vertex_prototypes[4], _vertex_prototypes[12]),
        Edge(18, _vertex_prototypes[7], _vertex_prototypes[12]),
        Edge(19, _vertex_prototypes[6], _vertex_prototypes[12]),

        Edge(20, _vertex_prototypes[1], _vertex_prototypes[13]),
        Edge(21, _vertex_prototypes[5], _vertex_prototypes[13]),
        Edge(22, _vertex_prototypes[6], _vertex_prototypes[13]),
        Edge(23, _vertex_prototypes[2], _vertex_prototypes[13]),
    )

    _num_windows = 4
    _num_window_types = 1
