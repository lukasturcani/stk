"""
Three Plus Six
==============

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class ThreePlusSix(Cage):
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
            topology_graph=stk.cage.ThreePlusSix((bb1, bb2)),
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
            topology_graph=stk.cage.ThreePlusSix(
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

        | 4-functional groups: 0 to 2
        | 2-functional groups: 3 to 8

    See :class:`.Cage` for more details and examples.

    """

    _x = 1
    _vertex_prototypes = (
        NonLinearVertex(0, [-2*_x, -_x*np.sqrt(3), 0], False),
        NonLinearVertex(1, [2*_x, -_x*np.sqrt(3), 0], False),
        NonLinearVertex(2, [0, _x*np.sqrt(3), 0], False),

        LinearVertex(3, [0, -2*_x*np.sqrt(3), _x], False),
        LinearVertex(4, [0, -2*_x*np.sqrt(3), -_x], False),

        LinearVertex(5, [2*_x, 0, _x], False),
        LinearVertex(6, [2*_x, 0, -_x], False),

        LinearVertex(7, [-2*_x, 0, _x], False),
        LinearVertex(8, [-2*_x, 0, -_x], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[3], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[3], _vertex_prototypes[1]),

        Edge(2, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[4], _vertex_prototypes[1]),

        Edge(4, _vertex_prototypes[5], _vertex_prototypes[1]),
        Edge(5, _vertex_prototypes[5], _vertex_prototypes[2]),

        Edge(6, _vertex_prototypes[6], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[6], _vertex_prototypes[2]),

        Edge(8, _vertex_prototypes[7], _vertex_prototypes[0]),
        Edge(9, _vertex_prototypes[7], _vertex_prototypes[2]),

        Edge(10, _vertex_prototypes[8], _vertex_prototypes[0]),
        Edge(11, _vertex_prototypes[8], _vertex_prototypes[2]),

    )

    _num_windows = 5
    _num_window_types = 2
