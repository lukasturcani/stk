"""
Four Plus Eight
===============

"""

from stk._internal.topology_graphs.cage.cage import Cage
from stk._internal.topology_graphs.cage.vertices import (
    LinearVertex,
    NonLinearVertex,
)
from stk._internal.topology_graphs.edge import Edge


class FourPlusEight(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='C1=CC=C(C(=C1)Br)Br',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='C1(=C(C(=C1Br)Br)Br)Br',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusEight((bb1, bb2)),
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
            smiles='C1=CC=C(C(=C1)Br)Br',
            functional_groups=[stk.BromoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='C1(=C(C(=C1Br)Br)Br)Br',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusEight(
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

        | 4-functional groups: 0 to 3
        | 2-functional groups: 4 to 11

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [-1, -1, 0], False),
        NonLinearVertex(1, [-1, 1, 0], False),
        NonLinearVertex(2, [1, -1, 0], False),
        NonLinearVertex(3, [1, 1, 0], False),
        LinearVertex(4, [-2, 0, 1], False),
        LinearVertex(5, [-2, 0, -1], False),
        LinearVertex(6, [0, 2, 1], False),
        LinearVertex(7, [0, 2, -1], False),
        LinearVertex(8, [0, -2, 1], False),
        LinearVertex(9, [0, -2, -1], False),
        LinearVertex(10, [2, 0, 1], False),
        LinearVertex(11, [2, 0, -1], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(4, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[3]),
        Edge(8, _vertex_prototypes[8], _vertex_prototypes[0]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[2]),
        Edge(10, _vertex_prototypes[9], _vertex_prototypes[0]),
        Edge(11, _vertex_prototypes[9], _vertex_prototypes[2]),
        Edge(12, _vertex_prototypes[10], _vertex_prototypes[2]),
        Edge(13, _vertex_prototypes[10], _vertex_prototypes[3]),
        Edge(14, _vertex_prototypes[11], _vertex_prototypes[2]),
        Edge(15, _vertex_prototypes[11], _vertex_prototypes[3]),
    )

    _num_windows = 6
    _num_window_types = 2
