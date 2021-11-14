"""
Two Plus Three
==============

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class TwoPlusThree(Cage):
    """
    Represents a cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='NC1CCCCC1N',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='O=Cc1cc(C=O)cc(C=O)c1',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.TwoPlusThree((bb1, bb2)),
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

    :class:`.MCHammer` optimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb1 = stk.BuildingBlock(
            smiles='NC1CCCCC1N',
            functional_groups=[stk.PrimaryAminoFactory()],
        )
        bb2 = stk.BuildingBlock(
            smiles='O=Cc1cc(C=O)cc(C=O)c1',
            functional_groups=[stk.AldehydeFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.TwoPlusThree(
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

        | 3-functional groups: (0, 1)
        | 2-functional groups: (2, 3, 4)

    See :class:`.Cage` for more details and examples.

    """

    _vertex_prototypes = (
        NonLinearVertex(0, [0, 0, 1]),
        NonLinearVertex(1, [0, 0, -1]),

        LinearVertex(2, [-1, -0.5*np.sqrt(3), 0], False),
        LinearVertex(3, [1, -0.5*np.sqrt(3), 0], False),
        LinearVertex(4, [0, 0.5*np.sqrt(3), 0], False),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(5, _vertex_prototypes[1], _vertex_prototypes[4]),
    )

    _num_windows = 3
    _num_window_types = 1
