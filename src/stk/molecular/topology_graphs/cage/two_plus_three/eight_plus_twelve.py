"""
Eight Plus Twelve
=================

"""

from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex
from ...topology_graph import Edge


class EightPlusTwelve(Cage):
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
            topology_graph=stk.cage.EightPlusTwelve((bb1, bb2)),
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

        | 3-functional groups: 0 to 7
        | 2-functional groups: 8 to 19

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [-1, 1, -1]),
        NonLinearVertex(1, [-1, -1, -1]),
        NonLinearVertex(2, [1, 1, -1]),
        NonLinearVertex(3, [1, -1, -1]),

        NonLinearVertex(4, [-1, 1, 1]),
        NonLinearVertex(5, [-1, -1, 1]),
        NonLinearVertex(6, [1, 1, 1]),
        NonLinearVertex(7, [1, -1, 1])
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        LinearVertex.init_at_center(
            id=8,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[2]),
        ),
        LinearVertex.init_at_center(
            id=9,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        LinearVertex.init_at_center(
            id=10,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[3]),
        ),
        LinearVertex.init_at_center(
            id=11,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[3]),
        ),

        LinearVertex.init_at_center(
            id=12,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[6]),
        ),
        LinearVertex.init_at_center(
            id=13,
            vertices=(_vertex_prototypes[4], _vertex_prototypes[5]),
        ),
        LinearVertex.init_at_center(
            id=14,
            vertices=(_vertex_prototypes[5], _vertex_prototypes[7]),
        ),
        LinearVertex.init_at_center(
            id=15,
            vertices=(_vertex_prototypes[6], _vertex_prototypes[7]),
        ),

        LinearVertex.init_at_center(
            id=16,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[4]),
        ),
        LinearVertex.init_at_center(
            id=17,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[5]),
        ),
        LinearVertex.init_at_center(
            id=18,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[6]),
        ),
        LinearVertex.init_at_center(
            id=19,
            vertices=(_vertex_prototypes[3], _vertex_prototypes[7]),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[8], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[8], _vertex_prototypes[2]),

        Edge(2, _vertex_prototypes[9], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[9], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[10], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[10], _vertex_prototypes[3]),

        Edge(6, _vertex_prototypes[11], _vertex_prototypes[2]),
        Edge(7, _vertex_prototypes[11], _vertex_prototypes[3]),

        Edge(8, _vertex_prototypes[12], _vertex_prototypes[4]),
        Edge(9, _vertex_prototypes[12], _vertex_prototypes[6]),

        Edge(10, _vertex_prototypes[13], _vertex_prototypes[4]),
        Edge(11, _vertex_prototypes[13], _vertex_prototypes[5]),

        Edge(12, _vertex_prototypes[14], _vertex_prototypes[5]),
        Edge(13, _vertex_prototypes[14], _vertex_prototypes[7]),

        Edge(14, _vertex_prototypes[15], _vertex_prototypes[6]),
        Edge(15, _vertex_prototypes[15], _vertex_prototypes[7]),

        Edge(16, _vertex_prototypes[16], _vertex_prototypes[0]),
        Edge(17, _vertex_prototypes[16], _vertex_prototypes[4]),

        Edge(18, _vertex_prototypes[17], _vertex_prototypes[1]),
        Edge(19, _vertex_prototypes[17], _vertex_prototypes[5]),

        Edge(20, _vertex_prototypes[18], _vertex_prototypes[2]),
        Edge(21, _vertex_prototypes[18], _vertex_prototypes[6]),

        Edge(22, _vertex_prototypes[19], _vertex_prototypes[3]),
        Edge(23, _vertex_prototypes[19], _vertex_prototypes[7]),

    )

    _num_windows = 6
    _num_window_types = 1
