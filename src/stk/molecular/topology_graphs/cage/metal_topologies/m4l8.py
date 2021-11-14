"""
M4L8
====

"""

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class M4L8(Cage):
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
            topology_graph=stk.cage.M4L8(
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
            topology_graph=stk.cage.M4L8(
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

    Metal building blocks with four functional groups are
    required for this topology.

    Ligand building blocks with two functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 4-functional groups: 0 to 3
        | 2-functional groups: 4 to 11

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [1, 0, 0]),
        NonLinearVertex(1, [0, 1, 0]),
        NonLinearVertex(2, [-1, 0, 0]),
        NonLinearVertex(3, [0, -1, 0]),

        LinearVertex(4, [1, 1, 0.5], False),
        LinearVertex(5, [1, 1, -0.5], False),

        LinearVertex(6, [1, -1, 0.5], False),
        LinearVertex(7, [1, -1, -0.5], False),

        LinearVertex(8, [-1, -1, 0.5], False),
        LinearVertex(9, [-1, -1, -0.5], False),

        LinearVertex(10, [-1, 1, 0.5], False),
        LinearVertex(11, [-1, 1, -0.5], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[5]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[6]),
        Edge(3, _vertex_prototypes[0], _vertex_prototypes[7]),

        Edge(4, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(6, _vertex_prototypes[1], _vertex_prototypes[10]),
        Edge(7, _vertex_prototypes[1], _vertex_prototypes[11]),

        Edge(8, _vertex_prototypes[2], _vertex_prototypes[10]),
        Edge(9, _vertex_prototypes[2], _vertex_prototypes[11]),
        Edge(10, _vertex_prototypes[2], _vertex_prototypes[8]),
        Edge(11, _vertex_prototypes[2], _vertex_prototypes[9]),

        Edge(12, _vertex_prototypes[3], _vertex_prototypes[8]),
        Edge(13, _vertex_prototypes[3], _vertex_prototypes[9]),
        Edge(14, _vertex_prototypes[3], _vertex_prototypes[6]),
        Edge(15, _vertex_prototypes[3], _vertex_prototypes[7]),
    )

    _num_windows = 2
    _num_window_types = 1
