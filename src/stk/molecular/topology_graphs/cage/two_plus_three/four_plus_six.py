"""
Four Plus Six
=============

"""

import numpy as np

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import LinearVertex, NonLinearVertex


class FourPlusSix(Cage):
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
            topology_graph=stk.cage.FourPlusSix((bb1, bb2)),
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
            topology_graph=stk.cage.FourPlusSix(
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

        | 3-functional groups: 0 to 3
        | 2-functional groups: 4 to 9

    See :class:`.Cage` for more details and examples.

    """

    # Vertices of a tetrahdron so that origin is at the origin. Source:
    # http://tinyurl.com/lc262h8.
    _v0, _v1, _v2, _v3 = _non_linears = (
        NonLinearVertex(0, [0, 0, np.sqrt(6)/2]),
        NonLinearVertex(1, [-1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        NonLinearVertex(2, [1, -np.sqrt(3)/3, -np.sqrt(6)/6]),
        NonLinearVertex(3, [0, 2*np.sqrt(3)/3, -np.sqrt(6)/6]),
    )

    _vertex_prototypes = (
        *_non_linears,
        LinearVertex.init_at_center(4, (_v0, _v1)),
        LinearVertex.init_at_center(5, (_v0, _v2)),
        LinearVertex.init_at_center(6, (_v0, _v3)),
        LinearVertex.init_at_center(7, (_v1, _v2)),
        LinearVertex.init_at_center(8, (_v1, _v3)),
        LinearVertex.init_at_center(9, (_v2, _v3)),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[4], _vertex_prototypes[0]),
        Edge(1, _vertex_prototypes[4], _vertex_prototypes[1]),
        Edge(2, _vertex_prototypes[5], _vertex_prototypes[0]),
        Edge(3, _vertex_prototypes[5], _vertex_prototypes[2]),
        Edge(4, _vertex_prototypes[6], _vertex_prototypes[0]),
        Edge(5, _vertex_prototypes[6], _vertex_prototypes[3]),
        Edge(6, _vertex_prototypes[7], _vertex_prototypes[1]),
        Edge(7, _vertex_prototypes[7], _vertex_prototypes[2]),
        Edge(8, _vertex_prototypes[8], _vertex_prototypes[1]),
        Edge(9, _vertex_prototypes[8], _vertex_prototypes[3]),
        Edge(10, _vertex_prototypes[9], _vertex_prototypes[2]),
        Edge(11, _vertex_prototypes[9], _vertex_prototypes[3])
    )

    _num_windows = 4
    _num_window_types = 1
