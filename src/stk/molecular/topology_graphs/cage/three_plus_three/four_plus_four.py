"""
Four Plus Four
==============

"""

from ...topology_graph import Edge
from ..cage import Cage
from ..vertices import NonLinearVertex


class FourPlusFour(Cage):
    """
    Represents a cube cage topology graph.

    Unoptimized construction

    .. moldoc::

        import moldoc.molecule as molecule
        import stk

        bb = stk.BuildingBlock(
            smiles='Brc1cc(Br)cc(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusFour((bb, )),
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

        bb = stk.BuildingBlock(
            smiles='Brc1cc(Br)cc(Br)c1',
            functional_groups=[stk.BromoFactory()],
        )
        cage = stk.ConstructedMolecule(
            topology_graph=stk.cage.FourPlusFour(
                building_blocks=(bb, ),
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

    Building blocks with three functional groups are required for
    this topology.

    When using a :class:`dict` for the `building_blocks` parameter,
    as in :ref:`cage-topology-graph-examples`:
    *Multi-Building Block Cage Construction*, a
    :class:`.BuildingBlock`, with the following number of functional
    groups, needs to be assigned to each of the following vertex ids:

        | 3-functional groups: 0 to 7

    See :class:`.Cage` for more details and examples.

    """

    _x = 1
    _vertex_prototypes = (
        NonLinearVertex(0, [-_x, _x, -_x], False),
        NonLinearVertex(1, [-_x, -_x, -_x], False),
        NonLinearVertex(2, [_x, _x, -_x], False),
        NonLinearVertex(3, [_x, -_x, -_x], False),

        NonLinearVertex(4, [-_x, _x, _x], False),
        NonLinearVertex(5, [-_x, -_x, _x], False),
        NonLinearVertex(6, [_x, _x, _x], False),
        NonLinearVertex(7, [_x, -_x, _x], False)
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[1]),
        Edge(1, _vertex_prototypes[0], _vertex_prototypes[2]),
        Edge(2, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[1], _vertex_prototypes[3]),
        Edge(4, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(6, _vertex_prototypes[2], _vertex_prototypes[3]),
        Edge(7, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(8, _vertex_prototypes[4], _vertex_prototypes[6]),
        Edge(9, _vertex_prototypes[4], _vertex_prototypes[5]),
        Edge(10, _vertex_prototypes[5], _vertex_prototypes[7]),
        Edge(11, _vertex_prototypes[6], _vertex_prototypes[7])
    )

    _num_windows = 6
    _num_window_types = 1
